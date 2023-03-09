using Pkg;
Pkg.activate(joinpath(@__DIR__, ".."));
using JLD2
using PGFPlotsX
using LaTeXStrings
using Colors
using LinearAlgebra

# Save plots as png or tikz
SAVEAS_PNG = false

color_list = distinguishable_colors(6, [RGB(1, 1, 1), RGB(0, 0, 0)], dropseed=true)
const def_linewidth = "ultra thick"
const color_KS = color_list[1]
const color_CW = color_list[2]
const color_YA = color_list[4]
const color_KGD = color_list[3]
const color_LIN = color_list[5]
# lineopts = @pgf {no_marks, "very thick", "error bars/y dir=both", "error bars/y explicit"}
lineopts = @pgf {no_marks, "ultra thick"}

data = load("data/rms_trajectory_error.jld2");
ranges = data["ranges"]
scales = data["scales"]

Δi_range = ranges.Δi_range
Δe_range = ranges.Δe_range
ΔM_range = ranges.ΔM_range
Δa_range = ranges.Δa_range

d_scale = scales.d_scale
t_scale = scales.t_scale

pos_err_M_cw = data["pos_err_M_cw"]
pos_err_M_ya = data["pos_err_M_ya"]
pos_err_M_kgd = data["pos_err_M_kgd"]
pos_err_M_ks = data["pos_err_M_ks"]
pos_err_M_lin = data["pos_err_M_lin"]
pos_err_i_cw = data["pos_err_i_cw"]
pos_err_i_ya = data["pos_err_i_ya"]
pos_err_i_kgd = data["pos_err_i_kgd"]
pos_err_i_ks = data["pos_err_i_ks"]
pos_err_i_lin = data["pos_err_i_lin"]
pos_err_e_cw = data["pos_err_e_cw"]
pos_err_e_ya = data["pos_err_e_ya"]
pos_err_e_kgd = data["pos_err_e_kgd"]
pos_err_e_ks = data["pos_err_e_ks"]
pos_err_e_lin = data["pos_err_e_lin"]
pos_err_a_cw = data["pos_err_a_cw"]
pos_err_a_ya = data["pos_err_a_ya"]
pos_err_a_kgd = data["pos_err_a_kgd"]
pos_err_a_ks = data["pos_err_a_ks"]
pos_err_a_lin = data["pos_err_a_lin"]

# scale errors to be in km
pos_err_M_cw = (d_scale / 1000) * pos_err_M_cw
pos_err_M_ya = (d_scale / 1000) * pos_err_M_ya
pos_err_M_kgd = (d_scale / 1000) * pos_err_M_kgd
pos_err_M_ks = (d_scale / 1000) * pos_err_M_ks
pos_err_M_lin = (d_scale / 1000) * pos_err_M_lin
pos_err_i_cw = (d_scale / 1000) * pos_err_i_cw
pos_err_i_ya = (d_scale / 1000) * pos_err_i_ya
pos_err_i_kgd = (d_scale / 1000) * pos_err_i_kgd
pos_err_i_ks = (d_scale / 1000) * pos_err_i_ks
pos_err_i_lin = (d_scale / 1000) * pos_err_i_lin
pos_err_e_cw = (d_scale / 1000) * pos_err_e_cw
pos_err_e_ya = (d_scale / 1000) * pos_err_e_ya
pos_err_e_kgd = (d_scale / 1000) * pos_err_e_kgd
pos_err_e_ks = (d_scale / 1000) * pos_err_e_ks
pos_err_e_lin = (d_scale / 1000) * pos_err_e_lin
pos_err_a_cw = (d_scale / 1000) * pos_err_a_cw
pos_err_a_ya = (d_scale / 1000) * pos_err_a_ya
pos_err_a_kgd = (d_scale / 1000) * pos_err_a_kgd
pos_err_a_ks = (d_scale / 1000) * pos_err_a_ks
pos_err_a_lin = (d_scale / 1000) * pos_err_a_lin

p_M = @pgf Axis(
    {
        xmajorgrids,
        ymajorgrids,
        xlabel = L"$\Delta$ Mean Anomaly (deg)",
        ylabel = "RMS Position error (km)",
        legend_pos = "north west",
        ymode = "log",
    },
    PlotInc({lineopts..., color = color_CW}, Coordinates(rad2deg.(ΔM_range), pos_err_M_cw)),
    PlotInc({lineopts..., color = color_YA, style = "dashed"}, Coordinates(rad2deg.(ΔM_range), pos_err_M_ya)),
    PlotInc({lineopts..., color = color_LIN, style = "dotted"}, Coordinates(rad2deg.(ΔM_range), pos_err_M_lin)),
    PlotInc({lineopts..., color = color_KGD, style = "dashed"}, Coordinates(rad2deg.(ΔM_range), pos_err_M_kgd)),
    PlotInc({lineopts..., color = color_KS}, Coordinates(rad2deg.(ΔM_range), pos_err_M_ks)),
)


p_i = @pgf Axis(
    {
        xmajorgrids,
        ymajorgrids,
        xlabel = L"$\Delta$ Inclination (deg)",
        ylabel = "RMS Position error (km)",
        legend_pos = "north west",
        ymode = "log"
    },
    PlotInc({lineopts..., color = color_CW}, Coordinates(rad2deg.(Δi_range), pos_err_i_cw)),
    PlotInc({lineopts..., color = color_YA, style = "dashed"}, Coordinates(rad2deg.(Δi_range), pos_err_i_ya)),
    PlotInc({lineopts..., color = color_LIN, style = "dotted"}, Coordinates(rad2deg.(Δi_range), pos_err_i_lin)),
    PlotInc({lineopts..., color = color_KGD, style = "dashed"}, Coordinates(rad2deg.(Δi_range), pos_err_i_kgd)),
    PlotInc({lineopts..., color = color_KS}, Coordinates(rad2deg.(Δi_range), pos_err_i_ks)),
)

p_e = @pgf LogLogAxis(
    {
        xmajorgrids,
        ymajorgrids,
        xlabel = L"$\Delta$ Eccentricity",
        ylabel = "RMS Position error (km)",
        legend_pos = "north west",
    },
    PlotInc({lineopts..., color = color_CW}, Coordinates(Δe_range, pos_err_e_cw)),
    PlotInc({lineopts..., color = color_YA, style = "dashed"}, Coordinates(Δe_range, pos_err_e_ya)),
    PlotInc({lineopts..., color = color_LIN, style = "dotted"}, Coordinates(Δe_range, pos_err_e_lin)),
    PlotInc({lineopts..., color = color_KGD, style = "dashed"}, Coordinates(Δe_range, pos_err_e_kgd)),
    PlotInc({lineopts..., color = color_KS}, Coordinates(Δe_range, pos_err_e_ks)),
)

p_a = @pgf Axis(
    {
        xmajorgrids,
        ymajorgrids,
        xlabel = L"$\Delta$ Semi-Major Axis (km)",
        ylabel = "RMS Position error (km)",
        legend_pos = "south east",
        ymode = "log"
    },
    PlotInc({lineopts..., color = color_CW}, Coordinates(Δa_range / 1000, pos_err_a_cw)),
    PlotInc({lineopts..., color = color_YA, style = "dashed"}, Coordinates(Δa_range / 1000, pos_err_a_ya)),
    PlotInc({lineopts..., color = color_LIN, style = "dotted"}, Coordinates(Δa_range / 1000, pos_err_a_lin)),
    PlotInc({lineopts..., color = color_KGD, style = "dashed"}, Coordinates(Δa_range / 1000, pos_err_a_kgd)),
    PlotInc({lineopts..., color = color_KS}, Coordinates(Δa_range / 1000, pos_err_a_ks)),
    SAVEAS_PNG ? Legend() : Legend(["CW", "YA", "LIN", "KGD", "KS"])
)

@pgf groupopts = {
    group_style = {group_size = "2 by 2", horizontal_sep = "0.8in", vertical_sep = "1.5cm"},
    height = "2in",
    width = "3.5in",
}
@pgf gp = GroupPlot(groupopts, p_M, p_i, p_e, p_a)


if SAVEAS_PNG
    pgfsave(joinpath("figs", "png", "rms_trajectory_error_mean_anomaly.png"), p_M, dpi=300)
    pgfsave(joinpath("figs", "png", "rms_trajectory_error_inclination.png"), p_i, dpi=300)
    pgfsave(joinpath("figs", "png", "rms_trajectory_error_eccentricity.png"), p_e, dpi=300)
    pgfsave(joinpath("figs", "png", "rms_trajectory_error_sma.png"), p_a, dpi=300)
    pgfsave(joinpath("figs", "png", "rms_trajectory_error.png"), gp, dpi=300)
else
    pgfsave(joinpath("figs", "rms_trajectory_error.tikz"), gp, include_preamble=false)
    pgfsave(joinpath("figs", "png", "rms_trajectory_error.png"), gp, dpi=300)
end

