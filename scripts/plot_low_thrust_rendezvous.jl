using Pkg;
Pkg.activate(joinpath(@__DIR__, ".."));
using JLD2
using PGFPlotsX
using LaTeXStrings
using Colors
using LinearAlgebra

# Save plots as png or tikz
SAVEAS_PDF = true

color_list = distinguishable_colors(6, [RGB(1, 1, 1), RGB(0, 0, 0)], dropseed=true)
const color_main = color_list[2]
const color_uR = color_main
const color_uT = color_list[3]
const color_uN = color_list[4]

# lineopts = @pgf {no_marks, "very thick", "error bars/y dir=both", "error bars/y explicit"}
lineopts = @pgf {no_marks, "very thick"}

data = load("data/low_thrust_rendezvous.jld2");
scales = data["scales"]

d_scale = scales.d_scale
t_scale = scales.t_scale

scales = data["scales"]
controls_eci_N = data["controls_eci_N"]
target_state_eci = data["target_state_eci"]
chaser_state_eci = data["chaser_state_eci"]
relative_states_rtn = data["relative_states_rtn"]
times = data["times"]

Ndata = length(times)

position_errors_eci_km = (d_scale / 1000) * [norm(chaser_state_eci[k][1:3] - target_state_eci[k][1:3]) for k = 1:Ndata]
velocity_errors_eci_km = (d_scale / (t_scale)) * [norm(chaser_state_eci[k][4:6] - target_state_eci[k][4:6]) for k = 1:Ndata]

controls_X_uN = 1e6 .* [controls_eci_N[k][1] for k = 1:Ndata]
controls_Y_uN = 1e6 .* [controls_eci_N[k][2] for k = 1:Ndata]
controls_Z_uN = 1e6 .* [controls_eci_N[k][3] for k = 1:Ndata]

relative_position_R_km = (d_scale / 1000) * [relative_states_rtn[k][1] for k = 1:Ndata]
relative_position_T_km = (d_scale / 1000) * [relative_states_rtn[k][2] for k = 1:Ndata]
relative_position_N_km = (d_scale / 1000) * [relative_states_rtn[k][3] for k = 1:Ndata]

p_err = @pgf Axis(
    {
        xmajorgrids,
        ymajorgrids,
        xlabel = "Time (orbits)",
        ylabel = "Error",
        legend_pos = "north east",
    },
    PlotInc({lineopts..., color = color_main}, Coordinates(times, position_errors_eci_km)),
    PlotInc({lineopts..., color = color_uN}, Coordinates(times, velocity_errors_eci_km)),
    Legend(["Position Error (km)", "Velocity Error (m/s)"])
)


p_controls = @pgf Axis(
    {
        xmajorgrids,
        ymajorgrids,
        xlabel = "Time (orbits)",
        ylabel = L"Thrust ($\mu \mathrm{m}/\mathrm{s}^2$)",
        legend_pos = "north east",
    },
    PlotInc({lineopts..., color = color_uR}, Coordinates(times, controls_X_uN)),
    PlotInc({lineopts..., color = color_uT}, Coordinates(times, controls_Y_uN)),
    PlotInc({lineopts..., color = color_uN}, Coordinates(times, controls_Z_uN)),
    Legend([L"\delta u_R", L"\delta u_T", L"\delta u_N"])
)

p_RT = @pgf Axis(
    {
        xmajorgrids,
        ymajorgrids,
        xlabel = "Tangential (km)",
        ylabel = "Radial (km)",
    },
    PlotInc({lineopts..., color = color_main}, Coordinates(relative_position_T_km, relative_position_R_km)),
)

p_RN = @pgf Axis(
    {
        xmajorgrids,
        ymajorgrids,
        xlabel = "Normal (km)",
        ylabel = "Radial (km)",
    },
    PlotInc({lineopts..., color = color_main}, Coordinates(relative_position_N_km, relative_position_R_km)),
)

@pgf groupopts = {
    group_style = {group_size = "2 by 2", horizontal_sep = "0.8in", vertical_sep = "1.5cm"},
    height = "2in",
    width = "3.5in",
}
@pgf gp = GroupPlot(groupopts, p_err, p_controls, p_RT, p_RN)

if SAVEAS_PDF
    pgfsave(joinpath("figs", "pdf", "low_thrust_rendezvous_position_error.pdf"), p_err, dpi=300)
    pgfsave(joinpath("figs", "pdf", "low_thrust_rendezvous_controls.pdf"), p_controls, dpi=300)
    pgfsave(joinpath("figs", "pdf", "low_thrust_rendezvous_radial_tangential.pdf"), p_RT, dpi=300)
    pgfsave(joinpath("figs", "pdf", "low_thrust_rendezvous_radial_normal.pdf"), p_RN, dpi=300)
    pgfsave(joinpath("figs", "pdf", "low_thrust_rendezvous_quad.pdf"), gp, dpi=300)
else
    pgfsave(joinpath("figs", "low_thrust_rendezvous.tikz"), gp, include_preamble=false)
end
