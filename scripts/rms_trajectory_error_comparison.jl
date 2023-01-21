using Pkg;
Pkg.activate(joinpath(@__DIR__, ".."));

using LinearAlgebra
using SatelliteDynamics
using Plots
using ProgressBars

# plot configuration
gr();
theme(:bright, linewidth=3, linealpha=0.9, legend_font_pointsize=10, palette=:seaborn_bright)
savefigures = false
savedata = true

SD = SatelliteDynamics

include("../src/state_transition_matrices.jl")
include("../src/propagate_trajectories.jl")


""" `osc_to_scaled(x_osc, d_scale, t_scale)`
convert from osc to scaled cartesian

    `x_osc = [sma, e, i, Ω, ω, M]`
"""
function osc_to_scaled(x_osc, d_scale, t_scale)

    # convert to cartesian
    x_cart = SD.sOSCtoCART(x_osc, use_degrees=false)
    r = x_cart[1:3] / d_scale
    v = x_cart[4:6] / (d_scale / t_scale)
    x_scaled = [r; v]

    replace!(x_scaled, -0.0 => 0.0)
    return x_scaled
end

# chief orbit (500km equatorial)
sma_c = 750e3 + SD.R_EARTH # semi-major axis
e_c = 0.0 # eccentricity
i_c = deg2rad(98.2) # inclination
ω_c = deg2rad(0.0) # argument of periapsis
Ω_c = deg2rad(0.0) # right angle of the ascending node
M_c = deg2rad(0.0) # Mean Anomaly

# base orbit period
T = 2 * pi * sqrt(sma_c^3 / SD.GM_EARTH)

# scaling
d_scale = sma_c
t_scale = T

# convert to scaled cartesian
osc0_chief = [sma_c, e_c, i_c, Ω_c, ω_c, M_c]
x0_chief_scaled = osc_to_scaled(osc0_chief, d_scale, t_scale)

GM_scaled = SD.GM_EARTH / (d_scale^3 / t_scale^2)
h = GM_scaled / norm(x0_chief_scaled[1:3]) - 0.5 * norm(x0_chief_scaled[4:6])^2

sma_scaled = sma_c / d_scale

# Set up j2 dynamics
sat_model_j2 = SatelliteModel( # No Drag, With J2
    1.0, # satellite_mass::Float64, -- Unused, but should be nonzero
    0.0, # drag_area::Float64, -- No Drag
    d_scale, # distance_scale::Float64, 
    t_scale, # time_scale::Float64
    1.0, # u_scale::Float64
    1.0, # rel_scale
    # GM::Float64=SatelliteDynamics.GM_EARTH,
    # J2::Float64=SatelliteDynamics.J2_EARTH,
    # R_earth::Float64=SatelliteDynamics.R_EARTH,
    # rho::Float64=1e-12, # TODO: need better drag model
    # CD::Float64=1.0,
    add_perturbations=true
)

# number of steps to take across Δ*_ranges
num_Δ_steps = 101 # make odd to get 0.0 for symmetric ranges

# test ranges
Δi_bounds = [deg2rad(-5.0), deg2rad(5.0)]
Δe_log10_bounds = [-4.0, -0.2] # log10 of eccentricity
ΔM_bounds = [deg2rad(-10.0), deg2rad(10.0)]
Δa_bounds = [-20.0e3, 20.0e3]

# increase sampling density around 0
cubic_sample_range(min, max, len) = ((range(-1, 1, length=len) .^ 3) .* (max - min) ./ 2) .+ (max + min) / 2

Δi_range = cubic_sample_range(Δi_bounds[1], Δi_bounds[2], num_Δ_steps)
Δe_range = 10 .^ range(Δe_log10_bounds[1], Δe_log10_bounds[2], length=num_Δ_steps)
ΔM_range = cubic_sample_range(ΔM_bounds[1], ΔM_bounds[2], num_Δ_steps)
Δa_range = cubic_sample_range(Δa_bounds[1], Δa_bounds[2], num_Δ_steps)

# Set up perturbed dynamics
sat_model_j2 = SatelliteModel( # No Drag, With J2
    1.0, # satellite_mass::Float64, -- Unused, but should be nonzero
    0.0, # drag_area::Float64, -- No Drag
    d_scale, # distance_scale::Float64, 
    t_scale, # time_scale::Float64
    1.0, # u_scale::Float64
    1.0, # rel_scale
    # GM::Float64=SatelliteDynamics.GM_EARTH,
    # J2::Float64=SatelliteDynamics.J2_EARTH,
    # R_earth::Float64=SatelliteDynamics.R_EARTH,
    # rho::Float64=1e-12, # TODO: need better drag model
    # CD::Float64=1.0,
    add_perturbations=true
)

# times to sample trajectories at
sample_times = range(0, 1.0, length=201)

############################################################
## Test different M
println("Evaluating Variations in Mean Anomaly")
pos_err_M_cw = zeros(num_Δ_steps)
pos_err_M_ya = zeros(num_Δ_steps)
pos_err_M_kgd = zeros(num_Δ_steps)
pos_err_M_ks = zeros(num_Δ_steps)
pos_err_M_lin = zeros(num_Δ_steps)

for k = ProgressBar(1:num_Δ_steps)
    ΔM = ΔM_range[k]
    osc_deputy = [sma_c, e_c, i_c, Ω_c, ω_c, ΔM + M_c]
    x_deputy_eci = osc_to_scaled(osc_deputy, d_scale, t_scale)

    (x_rel_traj_eci_j2,
        x_rel_traj_eci_cw,
        x_rel_traj_eci_ya,
        x_rel_traj_eci_kgd,
        x_rel_traj_eci_ks,
        x_rel_traj_eci_lin) = propagate_relative_state_all(x_deputy_eci, x0_chief_scaled, sample_times, GM_scaled, sat_model_j2)
    pos_err_M_cw[k] = relative_state_rms_position_error(x_rel_traj_eci_j2, x_rel_traj_eci_cw)
    pos_err_M_ya[k] = relative_state_rms_position_error(x_rel_traj_eci_j2, x_rel_traj_eci_ya)
    pos_err_M_kgd[k] = relative_state_rms_position_error(x_rel_traj_eci_j2, x_rel_traj_eci_kgd)
    pos_err_M_ks[k] = relative_state_rms_position_error(x_rel_traj_eci_j2, x_rel_traj_eci_ks)
    pos_err_M_lin[k] = relative_state_rms_position_error(x_rel_traj_eci_j2, x_rel_traj_eci_lin)
end

begin
    p = plot(rad2deg.(ΔM_range), (d_scale / 1000) * pos_err_M_ks, label="KS")
    plot!(rad2deg.(ΔM_range), (d_scale / 1000) * pos_err_M_cw, label="CW")
    plot!(rad2deg.(ΔM_range), (d_scale / 1000) * pos_err_M_ya, label="YA", style=:dash)
    plot!(rad2deg.(ΔM_range), (d_scale / 1000) * pos_err_M_kgd, label="KGD", style=:dot)
    plot!(rad2deg.(ΔM_range), (d_scale / 1000) * pos_err_M_lin, label="LIN", style=:dashdot)
    plot!(xlabel="Δ Mean Anomaly (deg)", ylabel="RMS Position Error (km)")
    plot!(yscale=:log10, ylims=[1e-8, 1e4])
    if savefigures
        savefig(p, "figs/rms_trajectory_error_mean_anomaly.pdf")
    end
    return p
end


############################################################
## Test different i
println("Evaluating Variations in Inclination")
pos_err_i_cw = zeros(num_Δ_steps)
pos_err_i_ya = zeros(num_Δ_steps)
pos_err_i_kgd = zeros(num_Δ_steps)
pos_err_i_ks = zeros(num_Δ_steps)
pos_err_i_lin = zeros(num_Δ_steps)

for k = ProgressBar(1:num_Δ_steps)
    Δi = Δi_range[k]
    osc_deputy = [sma_c, e_c, Δi + i_c, Ω_c, ω_c, M_c]
    x_deputy_eci = osc_to_scaled(osc_deputy, d_scale, t_scale)

    x_deputy_eci = osc_to_scaled(osc_deputy, d_scale, t_scale)

    (x_rel_traj_eci_j2,
        x_rel_traj_eci_cw,
        x_rel_traj_eci_ya,
        x_rel_traj_eci_kgd,
        x_rel_traj_eci_ks,
        x_rel_traj_eci_lin) = propagate_relative_state_all(x_deputy_eci, x0_chief_scaled, sample_times, GM_scaled, sat_model_j2)
    pos_err_i_cw[k] = relative_state_rms_position_error(x_rel_traj_eci_j2, x_rel_traj_eci_cw)
    pos_err_i_ya[k] = relative_state_rms_position_error(x_rel_traj_eci_j2, x_rel_traj_eci_ya)
    pos_err_i_kgd[k] = relative_state_rms_position_error(x_rel_traj_eci_j2, x_rel_traj_eci_kgd)
    pos_err_i_ks[k] = relative_state_rms_position_error(x_rel_traj_eci_j2, x_rel_traj_eci_ks)
    pos_err_i_lin[k] = relative_state_rms_position_error(x_rel_traj_eci_j2, x_rel_traj_eci_lin)
end

begin
    p = plot(rad2deg.(Δi_range), (d_scale / 1000) * pos_err_i_ks, label="KS")
    plot!(rad2deg.(Δi_range), (d_scale / 1000) * pos_err_i_cw, label="CW")
    plot!(rad2deg.(Δi_range), (d_scale / 1000) * pos_err_i_ya, label="YA", style=:dash)
    plot!(rad2deg.(Δi_range), (d_scale / 1000) * pos_err_i_kgd, label="KGD", style=:dot)
    plot!(rad2deg.(Δi_range), (d_scale / 1000) * pos_err_i_lin, label="LIN", style=:dashdot)
    plot!(xlabel="Δ Inclination (deg)", ylabel="Position Error (km)")
    plot!(yscale=:log10, ylims=[1e-12, 1e3])
    if savefigures
        savefig(p, "figs/rms_trajectory_error_inclination.pdf")
    end
    return p
end

############################################################
## Test different e
println("Evaluating Variations in Eccentricity")
pos_err_e_cw = zeros(num_Δ_steps)
pos_err_e_ya = zeros(num_Δ_steps)
pos_err_e_kgd = zeros(num_Δ_steps)
pos_err_e_ks = zeros(num_Δ_steps)
pos_err_e_lin = zeros(num_Δ_steps)

for k = ProgressBar(1:num_Δ_steps)
    Δe = Δe_range[k]
    # since we are interested in how eccentricity assumptions in the linearization affect the relative error,
    # we initialize the target and chaser at the same eccentricity, with slightly different positions
    osc_deputy = [sma_c, Δe, i_c + deg2rad(1e-3), Ω_c, ω_c, M_c + deg2rad(1e-3)] # approximately 1km separation
    osc_chief_e = [sma_c, Δe, i_c, Ω_c, ω_c, M_c]
    x_deputy_eci = osc_to_scaled(osc_deputy, d_scale, t_scale)
    x_chief_eci_e = osc_to_scaled(osc_chief_e, d_scale, t_scale)

    (x_rel_traj_eci_j2,
        x_rel_traj_eci_cw,
        x_rel_traj_eci_ya,
        x_rel_traj_eci_kgd,
        x_rel_traj_eci_ks,
        x_rel_traj_eci_lin) =
        propagate_relative_state_all(x_deputy_eci, x_chief_eci_e, sample_times, GM_scaled, sat_model_j2)
    pos_err_e_cw[k] = relative_state_rms_position_error(x_rel_traj_eci_j2, x_rel_traj_eci_cw)
    pos_err_e_ya[k] = relative_state_rms_position_error(x_rel_traj_eci_j2, x_rel_traj_eci_ya)
    pos_err_e_kgd[k] = relative_state_rms_position_error(x_rel_traj_eci_j2, x_rel_traj_eci_kgd)
    pos_err_e_ks[k] = relative_state_rms_position_error(x_rel_traj_eci_j2, x_rel_traj_eci_ks)
    pos_err_e_lin[k] = relative_state_rms_position_error(x_rel_traj_eci_j2, x_rel_traj_eci_lin)
end

begin
    p = plot(Δe_range, (d_scale / 1000) * pos_err_e_ks, label="KS")
    plot!(Δe_range, (d_scale / 1000) * pos_err_e_cw, label="CW")
    plot!(Δe_range, (d_scale / 1000) * pos_err_e_ya, label="YA", style=:dash)
    plot!(Δe_range, (d_scale / 1000) * pos_err_e_kgd, label="KGD", style=:dot)
    plot!(Δe_range, (d_scale / 1000) * pos_err_e_lin, label="LIN", style=:dashdot)
    plot!(xlabel="Eccentricity", ylabel="Position Error (km)")
    plot!(xscale=:log10, yscale=:log10, ylims=[1e-12, 1e3])
    if savefigures
        savefig(p, "figs/rms_trajectory_error_eccentricity.pdf")
    end
    return p
end

############################################################
## Test different sma
println("Evaluating Variations in Semi-Major Axis")
pos_err_a_cw = zeros(num_Δ_steps)
pos_err_a_ya = zeros(num_Δ_steps)
pos_err_a_kgd = zeros(num_Δ_steps)
pos_err_a_ks = zeros(num_Δ_steps)
pos_err_a_lin = zeros(num_Δ_steps)

for k = ProgressBar(1:num_Δ_steps)
    Δa = Δa_range[k]
    osc_deputy = [Δa + sma_c, e_c, i_c, Ω_c, ω_c, M_c]
    x_deputy_eci = osc_to_scaled(osc_deputy, d_scale, t_scale)

    (x_rel_traj_eci_j2,
        x_rel_traj_eci_cw,
        x_rel_traj_eci_ya,
        x_rel_traj_eci_kgd,
        x_rel_traj_eci_ks,
        x_rel_traj_eci_lin) =
        propagate_relative_state_all(x_deputy_eci, x0_chief_scaled, sample_times, GM_scaled, sat_model_j2)
    pos_err_a_cw[k] = relative_state_rms_position_error(x_rel_traj_eci_j2, x_rel_traj_eci_cw)
    pos_err_a_ya[k] = relative_state_rms_position_error(x_rel_traj_eci_j2, x_rel_traj_eci_ya)
    pos_err_a_kgd[k] = relative_state_rms_position_error(x_rel_traj_eci_j2, x_rel_traj_eci_kgd)
    pos_err_a_ks[k] = relative_state_rms_position_error(x_rel_traj_eci_j2, x_rel_traj_eci_ks)
    pos_err_a_lin[k] = relative_state_rms_position_error(x_rel_traj_eci_j2, x_rel_traj_eci_lin)
end

begin
    p = plot((1 / 1000) * Δa_range, (d_scale / 1000) * pos_err_a_ks, label="KS")
    plot!((1 / 1000) * Δa_range, (d_scale / 1000) * pos_err_a_cw, label="CW")
    plot!((1 / 1000) * Δa_range, (d_scale / 1000) * pos_err_a_ya, label="YA", style=:dash)
    plot!((1 / 1000) * Δa_range, (d_scale / 1000) * pos_err_a_kgd, label="KGD", style=:dot)
    plot!((1 / 1000) * Δa_range, (d_scale / 1000) * pos_err_a_lin, label="LIN", style=:dashdot)
    plot!(xlabel="Δ SMA (km) ", ylabel="Position Error (km)")
    plot!(yscale=:log10, ylims=[1e-10, 1e3])
    if savefigures
        savefig(p, "figs/rms_trajectory_error_sma.pdf")
    end
    return p
end


if savedata

    println("Saving Data...")
    ranges = (
        Δi_range=Δi_range,
        Δe_range=Δe_range,
        ΔM_range=ΔM_range,
        Δa_range=Δa_range
    )
    scales = (
        d_scale=d_scale,
        t_scale=t_scale
    )
    using JLD2
    jldsave("data/rms_trajectory_error.jld2";
        ranges=ranges,
        scales=scales,
        # M
        pos_err_M_cw=pos_err_M_cw,
        pos_err_M_ya=pos_err_M_ya,
        pos_err_M_kgd=pos_err_M_kgd,
        pos_err_M_ks=pos_err_M_ks,
        pos_err_M_lin=pos_err_M_lin,
        # i
        pos_err_i_cw=pos_err_i_cw,
        pos_err_i_ya=pos_err_i_ya,
        pos_err_i_kgd=pos_err_i_kgd,
        pos_err_i_ks=pos_err_i_ks,
        pos_err_i_lin=pos_err_i_lin,
        # e
        pos_err_e_cw=pos_err_e_cw,
        pos_err_e_ya=pos_err_e_ya,
        pos_err_e_kgd=pos_err_e_kgd,
        pos_err_e_ks=pos_err_e_ks,
        pos_err_e_lin=pos_err_e_lin,
        # a
        pos_err_a_cw=pos_err_a_cw,
        pos_err_a_ya=pos_err_a_ya,
        pos_err_a_kgd=pos_err_a_kgd,
        pos_err_a_ks=pos_err_a_ks,
        pos_err_a_lin=pos_err_a_lin
    )

end
