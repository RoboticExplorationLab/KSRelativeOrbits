using Pkg;
Pkg.activate(joinpath(@__DIR__, ".."));
using LinearAlgebra
using ForwardDiff
using SatelliteDynamics
using Plots
using OSQP
using SparseArrays
using StaticArrays
using Infiltrator
using JLD2

SD = SatelliteDynamics
theme(:bright, linewidth=3, linealpha=0.9, legend_font_pointsize=10, palette=:seaborn_bright)

savedata = true

include("../src/ks_transform.jl")
include("../src/trajopt_qp.jl")
include("../src/ks_trajopt_perturbed.jl")
include("../src/sim_dynamics.jl")

# Desired orbit initial conditions - ISS
sma_des = 417e3 + SD.R_EARTH
e_des = 0.0003
i_des = deg2rad(51.64)
ω_des = deg2rad(300)
Ω_des = deg2rad(0)
M_des = deg2rad(0.0)

xdes_osc = [sma_des, e_des, i_des, Ω_des, ω_des, M_des]
xdes_cart = SD.sOSCtoCART(xdes_osc, use_degrees=false)

# Initial Conditions

sma_0 = 415e3 + SD.R_EARTH
e_0 = 0.0004
i_0 = deg2rad(51.65)
ω_0 = deg2rad(300)
Ω_0 = deg2rad(0)
M_0 = deg2rad(-0.1)

x0_osc = [sma_0, e_0, i_0, Ω_0, ω_0, M_0]
x0_cart = SD.sOSCtoCART(x0_osc, use_degrees=false)

@show x0_cart .- xdes_cart
@show SD.sECItoRTN(xdes_cart, x0_cart)

# scaling parameters
d_scale = norm(xdes_cart[1:3])
period_des = 2 * pi * sqrt((sma_des)^3 / SD.GM_EARTH) # period
t_scale = period_des
GM_scaled = SD.GM_EARTH / (d_scale^3 / t_scale^2)

rdes = xdes_cart[1:3] / d_scale
vdes = xdes_cart[4:6] / (d_scale / t_scale)
xdes_scaled = [rdes; vdes]
qdes = state_cart_to_ks_newton(xdes_scaled)
hdes = GM_scaled / norm(rdes) - 0.5 * norm(vdes)^2

r0 = x0_cart[1:3] / d_scale
v0 = x0_cart[4:6] / (d_scale / t_scale)
x0_scaled = [r0; v0]
q0 = state_cart_to_ks_newton(x0_scaled; p_near=qdes[1:4])
h0 = GM_scaled / norm(r0) - 0.5 * norm(v0)^2

# state and control dimensions
Nx = length(q0) + 1
Nu = 3

# time stuff
ds = 0.05 # time step in KS state
Nk = 2000 # number of time steps

period_des_ks = 2 * pi / sqrt(hdes / 2)

# cost matrices
Q = 2.0 * I(Nx) / Nk
R = 1.0 * I(Nu) / Nk
Qf = 100.0 * Q

u_scale = 1.0 # / umax

rel_scale = norm(q0[1:4] .- qdes[1:4])

# Set up perturbed dynamics
sat_model = SatelliteModel( # No Drag, With J2
    1.0, # satellite_mass::Float64, -- Unused, but should be nonzero
    0.0, # drag_area::Float64, -- No Drag
    d_scale, # distance_scale::Float64, 
    t_scale, # time_scale::Float64
    u_scale,
    rel_scale, # rel_scale
    # GM::Float64=SatelliteDynamics.GM_EARTH,
    # J2::Float64=SatelliteDynamics.J2_EARTH,
    # R_earth::Float64=SatelliteDynamics.R_EARTH,
    # rho::Float64=1e-12, # TODO: need better drag model
    # CD::Float64=1.0,
    add_perturbations=true
)

# low thrust constraint
thrust_max = 20e-6 # Newtons
u_newtons_to_scaled_accel = (1.0 / sat_model.satellite_mass) * (t_scale^2) / d_scale
umax = thrust_max * u_newtons_to_scaled_accel
# thrust_max = d_scale * umax * sat_model.satellite_mass / (t_scale^2)

Nd_u = Nu * (Nk - 1)

function setup_low_thrust_constraints(qp, model)
    umax_scaled = umax * model.u_scale / model.rel_scale
    setup_u_bound_constraints!(qp, -umax_scaled, umax_scaled)
end

qh_sol_traj, u_sol_traj, q_des_traj, tvec = solve_relative_maneuver_perturbed(q0, qdes, ds, GM_scaled, Nk, Q, R, Qf, sat_model; num_constraints=Nd_u, constraints_setup_fn=setup_low_thrust_constraints)

# plot solution
pos_err = [d_scale * (norm(position_ks_to_inertial(qh_sol_traj[k][1:4] + q_des_traj[k][1:4]) - position_ks_to_inertial(q_des_traj[k][1:4]))) for k = 1:Nk]
orbits_ks = tvec / period_des_ks
begin
    plot(orbits_ks, pos_err)
    plot!(xlabel="Orbits", ylabel="Position Error (m)", title="Position Error")
end

rel_states_rtn = hcat((d_scale / 1000) * [Vector(SD.sECItoRTN(state_ks_to_inertial(q_des_traj[k][1:8]), position_ks_to_inertial(qh_sol_traj[k][1:4] + q_des_traj[k][1:4]))) for k = 1:Nk]...)
begin
    plot(rel_states_rtn[2, :], rel_states_rtn[1, :], title="LVLH Relative Position")
    plot!(xlabel="Local-Horizontal (km)", ylabel="Local-Vertical (km)")
end

begin
    plot(rel_states_rtn[3, :], rel_states_rtn[1, :], title="LVLH Relative Position")
    plot!(xlabel="Cross-Track (km)", ylabel="Local-Vertical (km)")
end
# plot controls
begin
    plot(orbits_ks[1:end-1], hcat(u_sol_traj...)', label=["u_R" "u_T" "u_N"])
    plot!(title="Controls")
end

# plot energy
begin
    plot(orbits_ks, [qh_sol_traj[k][9] for k = 1:Nk], title="Δ Energy", label="h̃")
end


## Nonlinear sim
function scaled_perturbed_dynamics(x_k, u_k, t)
    cartesian_J2_drag_perturbed_dynamics(x_k, u_k, t, sat_model)
end
# dynamics set up for using ODE.jl
# z = [x_chf, x_dep, u_dep]
function ode_dynamics!(zdot, z, p, t)
    zdot[1:6] .= scaled_perturbed_dynamics(z[1:6], zeros(3), t)

    # convert controls from chief RTN to ECI
    R = Matrix(SatelliteDynamics.rRTNtoECI(z[1:6]))
    u_eci = Vector(R * z[13:15])

    zdot[7:12] .= scaled_perturbed_dynamics(z[7:12], u_eci, t)
    zdot[13:15] .= 0.0 # controls
end

# condition for controls
function controls_condtion(out, z, t, integrator)
    p = integrator.p
    recompute_controls_out = t - (p.controls_update_period + p.controls_updated_time[1]) # condition for when to recompute controls
    controls_step_out = integrator.p.controls_times .- t # condition for when to get next control
    out .= [recompute_controls_out; controls_step_out]
end

# affect for controls
function controls_affect!(integrator, idx)
    p = integrator.p
    if idx == 1 # recompute controls
        update_time = integrator.t
        p.controls_updated_time[1] = update_time
        println("Recomputing Controls at t = $update_time")
        x_re_solve = integrator.u[7:12]
        x_re_solve_des = integrator.u[1:6]
        q_re_solve = state_cart_to_ks_newton(x_re_solve)
        q_re_solve_des = state_cart_to_ks_newton(x_re_solve_des; p_near=q_re_solve[1:4])
        _, controls_list, _, controls_times = solve_relative_maneuver_perturbed(q_re_solve, q_re_solve_des, ds, GM_scaled, Nk, Q, R, Qf, sat_model; num_constraints=Nd_u, constraints_setup_fn=setup_low_thrust_constraints)
        p.controls_list .= controls_list
        p.controls_times .= controls_times .+ update_time
        integrator.u[13:15] .= p.controls_list[1]
    else
        integrator.u[13:15] .= p.controls_list[idx-1]
    end
end


p = (controls_list=copy(u_sol_traj), controls_times=copy(tvec), controls_updated_time=[0.0], controls_update_period=1.0)
controls_cb = VectorContinuousCallback(controls_condtion, controls_affect!, length(p.controls_times) + 1)

tspan = (0.0, tvec[end-1])
z0 = [xdes_scaled; x0_scaled; u_sol_traj[1]]
prob = ODEProblem(ode_dynamics!, z0, tspan, p)
sol = solve(prob, Tsit5(), callback=controls_cb; abstol=1e-12, reltol=1e-13, saveat=0.1)

Nend = length(sol.u)

t = sol.t

pos_err = hcat(d_scale * [norm(sol.u[k][7:9] - sol.u[k][1:3]) for k = 1:Nend]...)
begin
    plot(t[1:Nend] * t_scale / 3600, pos_err[1:Nend], title="Chaser Position Error")
end

rel_states_rtn = hcat((d_scale / 1000) * [Vector(SD.sECItoRTN(sol.u[k][1:6], sol.u[k][7:9])) for k = 1:Nend]...)
begin
    plot(rel_states_rtn[2, :], rel_states_rtn[1, :], title="LVLH Relative Position")
    plot!(xlabel="Local-Horizontal (km)", ylabel="Local-Vertical (km)")
end

begin
    plot(rel_states_rtn[3, :], rel_states_rtn[1, :], title="LVLH Relative Position")
    plot!(xlabel="Cross-Track (km)", ylabel="Local-Vertical (km)")
end

begin
    plot(t[1:Nend] * t_scale / 3600, hcat([sol.u[k][13:15] for k = 1:Nend]...)', title="Chaser Controls")
end


if savedata
    Ndata = length(sol.u)

    controls_eci_N = [sol.u[k][13:15] for k = 1:Ndata] ./ u_newtons_to_scaled_accel
    target_state_eci = [sol.u[k][1:6] for k = 1:Ndata]
    chaser_state_eci = [sol.u[k][7:12] for k = 1:Ndata]
    relative_states_rtn = [Vector(SD.sECItoRTN(target_state_eci[k], chaser_state_eci[k])) for k = 1:Ndata]
    times = sol.t

    println("Saving Data...")
    scales = (
        d_scale=d_scale,
        t_scale=t_scale
    )
    jldsave("data/low_thrust_rendezvous.jld2";
        scales=scales,
        controls_eci_N=controls_eci_N,
        target_state_eci=target_state_eci,
        chaser_state_eci=chaser_state_eci,
        relative_states_rtn=relative_states_rtn,
        times=times
    )

end

time_solver = false
using BenchmarkTools
if time_solver
    @benchmark solve_relative_maneuver_perturbed(q0, qdes, ds, GM_scaled, Nk, Q, R, Qf, sat_model; num_constraints=Nd_u, constraints_setup_fn=setup_low_thrust_constraints)
end