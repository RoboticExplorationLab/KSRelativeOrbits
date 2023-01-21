using Pkg;
Pkg.activate(joinpath(@__DIR__, ".."));

using LinearAlgebra
using SatelliteDynamics
using Plots
using ProgressBars

SD = SatelliteDynamics

include("../src/state_transition_matrices.jl")
include("orbit_utils.jl")

function propagate_relative_state_all(x0_deputy, x0_chief, times, GM, sat_model::SatelliteModel)
    Nk = length(times)

    # propagate nonlinear chief and deputy states
    x_traj_chief_eci_j2 = propagate_j2_state(x0_chief, times, sat_model)
    x_traj_deputy_eci_j2 = propagate_j2_state(x0_deputy, times, sat_model)

    x_rel_traj_eci_j2 = [x_traj_deputy_eci_j2[k] .- x_traj_chief_eci_j2[k] for k = 1:Nk]

    # propagate STM states
    x_rel_traj_eci_cw = propagate_cw_state(x0_deputy, x_traj_chief_eci_j2, times, GM)
    x_rel_traj_eci_ya = propagate_ya_state(x0_deputy, x_traj_chief_eci_j2, times, GM)
    x_rel_traj_eci_kgd = propagate_kgd_state(x0_deputy, x_traj_chief_eci_j2, times, sat_model)

    # propagate KS state
    _, x_rel_traj_eci_ks, _ = propagate_ks_state(x0_deputy, x0_chief, times, sat_model)

    # propagate linearized eci state
    _, x_rel_traj_eci_lin = propagate_lin_eci_state(x0_deputy, x0_chief, times, sat_model)

    return x_rel_traj_eci_j2, x_rel_traj_eci_cw, x_rel_traj_eci_ya, x_rel_traj_eci_kgd, x_rel_traj_eci_ks, x_rel_traj_eci_lin
end

function propagate_j2_state(x0, times, sat_model::SatelliteModel)

    function perturbed_dynamics!(ẋ, x, p, t)
        ẋ .= cartesian_J2_drag_perturbed_dynamics(x, zeros(3), t, sat_model)
    end

    prob = ODEProblem(perturbed_dynamics!, x0, (times[1], times[end]))

    sol = solve(prob, Tsit5(); abstol=1e-12, reltol=1e-13, saveat=times)

    x_traj_kep = [sol.u[k][1:6] for k = 1:length(sol.u)]

    return x_traj_kep
end

function propagate_kep_state(x0, times, GM)

    function scaled_keplarian_dynamics!(ẋ, x, p, t)
        ẋ .= keplarian_dynamics(x[1:6], zeros(3), t; GM=GM) # no control inputs
    end

    prob = ODEProblem(scaled_keplarian_dynamics!, x0, (times[1], times[end]))

    sol = solve(prob, Tsit5(); abstol=1e-12, reltol=1e-13, saveat=times)

    x_traj_kep = [sol.u[k][1:6] for k = 1:length(sol.u)]

    return x_traj_kep
end

""" `propagate_cw_state(x0_deputy, x_traj_chief, times, GM_scaled)`
Propagate the relative state using each state transition matrix formulation.
"""
function propagate_cw_state(x0_deputy_eci, x_traj_chief_eci, times, GM_scaled)
    Nk = length(times)

    x0_chief_eci = x_traj_chief_eci[1]
    osc_chief = SD.sCARTtoOSC(x0_chief_eci, GM=GM_scaled)
    sma_c, e_c, i_c, Ω_c, ω_c, M_c = osc_chief

    x_rel_0_eci = x0_deputy_eci .- x0_chief_eci

    x_rel_traj_eci_cw = [zeros(6) for k = 1:Nk]

    # initial condition
    x_rel_traj_eci_cw[1] .= x_rel_0_eci

    for k = 1:Nk-1
        xk_chief_eci = x_traj_chief_eci[k]
        xkp1_chief_eci = x_traj_chief_eci[k+1]
        dt = times[k+1] - times[k]

        # CW
        xk_rel_eci_cw = x_rel_traj_eci_cw[k] # get current cw relative state in ECI
        xk_rel_rtn_cw = SD.sECItoRTN(xk_chief_eci, xk_chief_eci .+ xk_rel_eci_cw) # transform from ECI to RTN at current chief state
        xkp1_rel_rtn_cw = transition_relative_state_via_cw(xk_rel_rtn_cw, GM_scaled, sma_c, dt) # use state transition matrix to find next state
        x_rel_traj_eci_cw[k+1] .= SD.sRTNtoECI(xkp1_chief_eci, xkp1_rel_rtn_cw) .- xkp1_chief_eci # transform back to ECI at next chief state
    end

    return x_rel_traj_eci_cw
end

""" `propagate_ya_state(x0_deputy_eci, x0_chief_eci, times, GM)`
Propagate the relative state using each state transition matrix formulation.
"""
function propagate_ya_state(x0_deputy_eci, x_traj_chief_eci, times, GM)
    Nk = length(times)

    x_rel_traj_rtn_ya = [zeros(6) for k = 1:Nk]

    # initial condition
    x0_chief_eci = x_traj_chief_eci[1]
    x0_rel_rtn_ya = SD.sECItoRTN(x0_chief_eci, x0_deputy_eci)
    x_rel_traj_rtn_ya[1] .= x0_rel_rtn_ya

    osc0_chief = SD.sCARTtoOSC(x0_chief_eci, GM=GM)
    sma0_c, e0_c, i0_c, Ω0_c, ω0_c, M0_c = osc0_chief

    E0 = E_from_M(M0_c, e0_c)
    f0 = theta_from_E(E0, e0_c)

    r0 = norm(x0_chief_eci[1:3])
    p0 = sma0_c * (1 - e0_c^2)
    h0 = sqrt(GM * p0)
    n0 = sqrt(GM / (sma0_c^3))

    f0_dot = h0 / (r0^2)
    r0_dot = p0 * e0_c * f0_dot * sin(f0) / (1 + e0_c * cos(f0))^2

    # carter normalization
    x, y, z, vx, vy, vz = x0_rel_rtn_ya
    x0_rel_carter_ya = [
        x / r0,
        (r0 * vx - x * r0_dot) / (f0_dot * r0^2),
        y / r0,
        (r0 * vy - y * r0_dot) / (f0_dot * r0^2),
        z / r0,
        (r0 * vz - z * r0_dot) / (f0_dot * r0^2)]

    ϕ_0_inv = ya_phi_inv0(e0_c, f0)

    K_intermediate = ϕ_0_inv * x0_rel_carter_ya

    # time
    t0 = times[1]

    for k = 1:Nk
        Tk = times[k] - t0

        Mk = M0_c + Tk * n0
        Ek = E_from_M(Mk, e0_c)
        fk = theta_from_E(Ek, e0_c)

        ϕ_k = ya_phi(e0_c, fk, GM, h0, Tk)

        xk_carter = ϕ_k * K_intermediate

        rk = p0 / (1 + e0_c * cos(fk))
        fk_dot = h0 / (rk^2)
        rk_dot = p0 * e0_c * fk_dot * sin(fk) / (1 + e0_c * cos(fk))^2

        xbar, xbar_p, ybar, ybar_p, zbar, zbar_p = xk_carter
        xk_rel_rtn = [
            xbar * rk,
            ybar * rk,
            zbar * rk,
            xbar_p * fk_dot * rk + xbar * rk_dot,
            ybar_p * fk_dot * rk + ybar * rk_dot,
            zbar_p * fk_dot * rk + zbar * rk_dot,
        ]

        x_rel_traj_rtn_ya[k] .= xk_rel_rtn
    end

    # convert to relative ECI state
    x_rel_traj_eci_ya = [SD.sRTNtoECI(x_traj_chief_eci[k], x_rel_traj_rtn_ya[k]) .- x_traj_chief_eci[k] for k = 1:Nk]

    return x_rel_traj_eci_ya
end

""" `propagate_kgd_state(x0_deputy, x_traj_chief, times, sat_model)`
Propagate the relative state using each state transition matrix formulation.
"""
function propagate_kgd_state(x0_deputy_eci, x_traj_chief_eci, times, sat_model)
    Nk = length(times)

    x0_chief_eci = x_traj_chief_eci[1]
    x_rel_0_eci = x0_deputy_eci .- x0_chief_eci

    x_rel_traj_eci_kgd = [zeros(6) for k = 1:Nk]

    # initial condition
    x_rel_traj_eci_kgd[1] .= x_rel_0_eci

    t0 = times[1]

    for k = 1:Nk
        xk_chief_eci = x_traj_chief_eci[k]
        T = times[k] - t0

        # KGD
        x_rel_traj_eci_kgd[k] .= transition_relative_state_via_kgd(x_rel_0_eci, x0_chief_eci, xk_chief_eci, T, sat_model)
    end

    return x_rel_traj_eci_kgd
end

function relative_state_rms_position_error(x_rel_traj_true, x_rel_traj_comp)
    Nk = length(x_rel_traj_true)
    accum = 0.0
    for k = 1:Nk
        diff = x_rel_traj_comp[k][1:3] .- x_rel_traj_true[k][1:3]
        accum += diff'diff
    end
    return sqrt(accum / Nk)
end

function propagate_ks_state(x0_deputy_eci, x0_chief_eci, times, sat_model::SatelliteModel)

    t_start = times[1]
    dt = times[end] - times[end-1]
    t_end = times[end] + dt # make sure to end after times[end] to get full trajectory

    # convert from eci to ks
    x0_chief_ks = state_cart_to_ks_newton(x0_chief_eci)
    x0_deputy_ks = state_cart_to_ks_newton(x0_deputy_eci; p_near=x0_chief_ks[1:4])
    x0_rel_ks = x0_deputy_ks .- x0_chief_ks

    # scaling
    d_scale = sat_model.distance_scale
    t_scale = sat_model.time_scale
    GM_scaled = SD.GM_EARTH / (d_scale^3 / t_scale^2)

    # energy
    h0_chief = ks_h_energy(x0_chief_ks[1:4], x0_chief_ks[5:8], GM_scaled)
    h0_deputy = ks_h_energy(x0_deputy_ks[1:4], x0_deputy_ks[5:8], GM_scaled)
    h0_rel = h0_deputy - h0_chief

    # dynamics with single vector input for ForwardDiff
    Nu = 3
    ks_full_dynamics_vec(xh) = ks_full_dynamics(sat_model, xh, zeros(Nu))

    # state + stm dynamics
    function ks_full_and_relative_dynamics!(z_ds, z, p, s)
        xh_chief = z[1:9]
        xh_rel = z[10:18]
        # nonlinear dynamics for chief
        xh_chief_ds = ks_full_dynamics_vec(xh_chief)

        # linearize around chief position
        ∂f∂xh = ForwardDiff.jacobian(ks_full_dynamics_vec, xh_chief)

        # linear relative dynamics
        xh_rel_ds = ∂f∂xh * xh_rel

        t_chief_ds = xh_chief[1:4]'xh_chief[1:4]

        xh_deputy = xh_chief .+ xh_rel
        t_deputy_ds = xh_deputy[1:4]'xh_deputy[1:4]

        z_ds .= [xh_chief_ds; xh_rel_ds; t_chief_ds; t_deputy_ds]
    end

    # initial condition
    z_0 = [x0_chief_ks; h0_chief; x0_rel_ks; h0_rel; t_start; t_start]
    Nz = length(z_0)

    # save at specific real time points
    Nsaves = length(times)
    xh_traj_chief_ks = [zeros(9) for k = 1:Nsaves]
    xh_traj_deputy_ks = [zeros(9) for k = 1:Nsaves]
    xh_traj_chief_ks[1] .= [x0_chief_ks; h0_chief]
    xh_traj_deputy_ks[1] .= [x0_deputy_ks; h0_deputy]

    t_chief_traj = zeros(Nsaves)
    t_dep_traj = zeros(Nsaves)
    t_chief_traj[1] = t_start
    t_dep_traj[1] = t_start

    function time_condition_separate!(out, z, s, integrator)
        t_chief = z[end-1]
        t_deputy = z[end]
        out .= [times .- t_chief; times .- t_deputy]
    end

    function time_affect_separate!(integrator, idx)
        if idx <= Nsaves
            # chief timepoint
            xh_traj_chief_ks[idx] .= integrator.u[1:9]
            t_chief_traj[idx] = integrator.u[end-1]
        else
            # deputy timepoint
            idx -= Nsaves
            xh_traj_deputy_ks[idx] .= integrator.u[10:18] .+ integrator.u[1:9]
            t_dep_traj[idx] = integrator.u[end]
        end
    end

    function time_condition_together!(out, z, s, integrator)
        t_chief = z[end-1]
        out .= times .- t_chief
    end

    function time_affect_together!(integrator, idx)
        # chief timepoint
        xh_traj_chief_ks[idx] .= integrator.u[1:9]
        t_chief_traj[idx] = integrator.u[end-1]
        # deputy timepoint
        xh_traj_deputy_ks[idx] .= integrator.u[10:18] .+ integrator.u[1:9]
        t_dep_traj[idx] = integrator.u[end]
    end

    # solver has issues when relative state is zero
    if norm(x0_rel_ks[1:4]) > 1e-14
        time_cb = VectorContinuousCallback(time_condition_separate!, time_affect_separate!, 2 * Nsaves)
    else
        time_cb = VectorContinuousCallback(time_condition_together!, time_affect_together!, Nsaves)
    end

    prob = ODEProblem(ks_full_and_relative_dynamics!, z_0, (t_start, t_end))
    sol = solve(prob, Tsit5(); abstol=1e-12, reltol=1e-13, callback=time_cb)

    # convert solution to ECI
    x_traj_chief_eci_ks = [state_ks_to_cart(xh_traj_chief_ks[k][1:8]) for k = 1:Nsaves]
    x_traj_deputy_eci_ks = [state_ks_to_cart(xh_traj_deputy_ks[k][1:8]) for k = 1:Nsaves]
    x_rel_traj_eci_ks = [x_traj_deputy_eci_ks[k][1:6] .- x_traj_chief_eci_ks[k][1:6] for k = 1:Nsaves]

    return x_traj_chief_eci_ks, x_rel_traj_eci_ks, t_chief_traj, t_dep_traj
end

function propagate_lin_eci_state(x0_deputy_eci, x0_chief_eci, times, sat_model::SatelliteModel)

    t_start = times[1]
    dt = times[end] - times[end-1]
    t_end = times[end] + dt # make sure to end after times[end] to get full trajectory

    # relative state
    x0_rel_eci = x0_deputy_eci .- x0_chief_eci

    # dynamics with single vector input for ForwardDiff
    Nu = 3
    eci_full_dynamics_vec(x, t) = cartesian_J2_drag_perturbed_dynamics(x, zeros(Nu), t, sat_model)

    # state + stm dynamics
    function eci_full_and_relative_dynamics!(z_dot, z, p, t)
        x_chief = z[1:6]
        x_rel = z[7:12]
        # nonlinear dynamics for chief
        x_chief_dot = eci_full_dynamics_vec(x_chief, t)

        # linearize around chief position
        ∂f∂x = ForwardDiff.jacobian(x_ -> eci_full_dynamics_vec(x_, t), x_chief)

        # linear relative dynamics
        x_rel_dot = ∂f∂x * x_rel

        z_dot .= [x_chief_dot; x_rel_dot]
    end

    # initial condition
    z_0 = [x0_chief_eci; x0_rel_eci]
    Nz = length(z_0)

    # save at specific real time points
    Nsaves = length(times)
    x_traj_chief_eci_lin = [zeros(6) for k = 1:Nsaves]
    x_traj_rel_eci_lin = [zeros(6) for k = 1:Nsaves]
    x_traj_chief_eci_lin[1] .= x0_chief_eci
    x_traj_rel_eci_lin[1] .= x0_rel_eci

    function time_condition!(out, z, t, integrator)
        out .= times .- t
    end

    function time_affect!(integrator, idx)
        x_traj_chief_eci_lin[idx] .= integrator.u[1:6]
        x_traj_rel_eci_lin[idx] .= integrator.u[7:12]
    end

    time_cb = VectorContinuousCallback(time_condition!, time_affect!, Nsaves)

    prob = ODEProblem(eci_full_and_relative_dynamics!, z_0, (t_start, t_end))
    sol = solve(prob, Tsit5(); abstol=1e-12, reltol=1e-13, callback=time_cb)

    return x_traj_chief_eci_lin, x_traj_rel_eci_lin
end