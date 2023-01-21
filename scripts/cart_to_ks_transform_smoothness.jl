using Pkg;
Pkg.activate(joinpath(@__DIR__, ".."));

using LinearAlgebra
using SatelliteDynamics
using Plots
using ProgressBars

# plot configuration
gr();
theme(:bright, linewidth=3, linealpha=0.9, legend_font_pointsize=10, palette=:seaborn_bright)
savedata = true

SD = SatelliteDynamics

include("../src/ks_transform.jl")
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

sma_c = 1e5
e_c = 0 # eccentricity
i_c = 0
ω_c = 0
Ω_c = 0
M_c = 0

# base orbit period
T = 2 * pi * sqrt(sma_c^3 / SD.GM_EARTH)

# scaling
d_scale = sma_c
t_scale = T

# convert to scaled cartesian
osc0 = [sma_c, e_c, i_c, Ω_c, ω_c, M_c]
x0 = osc_to_scaled(osc0, d_scale, t_scale)

GM_scaled = SD.GM_EARTH / (d_scale^3 / t_scale^2)

times = 0:0.01:1.0

x_traj_kep = propagate_kep_state(x0, times, GM_scaled)

# convert to ks
Nk = length(times)
x_traj_ks_fixed = [zeros(4) for k = 1:Nk]
x_traj_ks_newton = [zeros(4) for k = 1:Nk]
p_near = [1.0, 0.0, 0, 0]

for k = 1:Nk

    x_traj_ks_fixed[k] .= position_cart_to_ks_plus(x_traj_kep[k][1:3])
    x_traj_ks_newton[k] .= position_cart_to_ks_newton(x_traj_kep[k][1:3]; p_near=p_near)
    p_near = x_traj_ks_newton[k]
end

begin
    plot(times, hcat(x_traj_ks_fixed...)', label="fixed")#, color=:red)
    plot!(times, hcat(x_traj_ks_newton...)', label="Newton", style=:dash)#, color=:blue)
end

if savedata
    using JLD2
    jldsave("data/cart_to_ks_transform_smoothness.jld2";
        times=collect(times),
        x_traj_ks_fixed=x_traj_ks_fixed,
        x_traj_ks_newton=x_traj_ks_newton)
end