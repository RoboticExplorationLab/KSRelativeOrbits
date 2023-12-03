using Pkg;
Pkg.activate(joinpath(@__DIR__, ".."));
using JLD2
using PGFPlotsX
using LaTeXStrings
using Colors
using LinearAlgebra

color_list = distinguishable_colors(6, [RGB(1, 1, 1), RGB(0, 0, 0)], dropseed=true)
const color_none = RGBA(0, 0, 0, 0)
const def_linewidth = "ultra thick"
const color_KS = color_list[1]
const color_CW = color_list[2]
const color_YA = color_list[4]
const color_KGD = color_list[3]
const color_LIN = color_list[5]

color_mode = "_dark_mode"
# color_mode = "" # normal

if color_mode == "_dark_mode"
    const color_grid = RGBA(([148, 148, 148, 255] ./ 255)...)
    const color_text = RGBA(([205, 209, 209, 255] ./ 255)...)
    const color_axis = color_text
else
    const color_grid = RGBA(([191, 191, 191, 255] ./ 255)...)
    const color_text = RGBA(([0, 0, 0, 255] ./ 255)...)
    const color_axis = color_text
end

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


p_M(opacities) = @pgf Axis(
    {
        "grid style" = {"color" = color_grid},
        "label style" = {"color" = color_text},
        "tick label style" = {"color" = color_text},
        "axis line style" = {"color" = color_axis},
        label_style = raw"font=\LARGE",
        xmajorgrids,
        ymajorgrids,
        xlabel = L"$\Delta$ Mean Anomaly (deg)",
        ylabel = "RMS Position error (km)",
        legend_pos = "north west",
        ymode = "log",
    },
    PlotInc({lineopts..., opacity = opacities[1], color = color_CW}, Coordinates(rad2deg.(ΔM_range), pos_err_M_cw)),
    PlotInc({lineopts..., opacity = opacities[2], color = color_YA, style = "dashed"}, Coordinates(rad2deg.(ΔM_range), pos_err_M_ya)),
    PlotInc({lineopts..., opacity = opacities[3], color = color_KGD, style = "dashed"}, Coordinates(rad2deg.(ΔM_range), pos_err_M_kgd)),
    # PlotInc({lineopts..., opacity = opacities[4], color = color_LIN, style = "dotted"}, Coordinates(rad2deg.(ΔM_range), pos_err_M_lin)),
    PlotInc({lineopts..., opacity = opacities[5], color = color_KS}, Coordinates(rad2deg.(ΔM_range), pos_err_M_ks)),
)

p_M_cw = p_M([1, 0, 0, 0, 0])
p_M_ya = p_M([1, 1, 0, 0, 0])
p_M_kgd = p_M([1, 1, 1, 0, 0])
p_M_lin = p_M([1, 1, 1, 1, 0])
p_M_all = p_M([1, 1, 1, 1, 1])


p_i(opacities) = @pgf Axis(
    {
        "grid style" = {"color" = color_grid},
        "label style" = {"color" = color_text},
        "tick label style" = {"color" = color_text},
        "axis line style" = {"color" = color_axis},
        label_style = raw"font=\LARGE",
        xmajorgrids,
        ymajorgrids,
        xlabel = L"$\Delta$ Inclination (deg)",
        ylabel = "RMS Position error (km)",
        legend_pos = "north west",
        ymode = "log"
    },
    PlotInc({lineopts..., opacity = opacities[1], color = color_CW}, Coordinates(rad2deg.(Δi_range), pos_err_i_cw)),
    PlotInc({lineopts..., opacity = opacities[2], color = color_YA, style = "dashed"}, Coordinates(rad2deg.(Δi_range), pos_err_i_ya)),
    PlotInc({lineopts..., opacity = opacities[3], color = color_KGD, style = "dashed"}, Coordinates(rad2deg.(Δi_range), pos_err_i_kgd)),
    # PlotInc({lineopts..., opacity = opacities[4], color = color_LIN, style = "dotted"}, Coordinates(rad2deg.(Δi_range), pos_err_i_lin)),
    PlotInc({lineopts..., opacity = opacities[5], color = color_KS}, Coordinates(rad2deg.(Δi_range), pos_err_i_ks)),
)
p_i_cw = p_i([1, 0, 0, 0, 0])
p_i_ya = p_i([1, 1, 0, 0, 0])
p_i_kgd = p_i([1, 1, 1, 0, 0])
p_i_lin = p_i([1, 1, 1, 1, 0])
p_i_all = p_i([1, 1, 1, 1, 1])

p_e(opacities) = @pgf LogLogAxis(
    {
        "grid style" = {"color" = color_grid},
        "label style" = {"color" = color_text},
        "tick label style" = {"color" = color_text},
        "axis line style" = {"color" = color_axis},
        label_style = raw"font=\LARGE",
        xmajorgrids,
        ymajorgrids,
        xlabel = L"$\Delta$ Eccentricity",
        ylabel = "RMS Position error (km)",
        legend_pos = "north west",
    },
    PlotInc({lineopts..., opacity = opacities[1], color = color_CW}, Coordinates(Δe_range, pos_err_e_cw)),
    PlotInc({lineopts..., opacity = opacities[2], color = color_YA, style = "dashed"}, Coordinates(Δe_range, pos_err_e_ya)),
    PlotInc({lineopts..., opacity = opacities[3], color = color_KGD, style = "dashed"}, Coordinates(Δe_range, pos_err_e_kgd)),
    # PlotInc({lineopts..., opacity = opacities[4], color = color_LIN, style = "dotted"}, Coordinates(Δe_range, pos_err_e_lin)),
    PlotInc({lineopts..., opacity = opacities[5], color = color_KS}, Coordinates(Δe_range, pos_err_e_ks)),
)
p_e_cw = p_e([1, 0, 0, 0, 0])
p_e_ya = p_e([1, 1, 0, 0, 0])
p_e_kgd = p_e([1, 1, 1, 0, 0])
p_e_lin = p_e([1, 1, 1, 1, 0])
p_e_all = p_e([1, 1, 1, 1, 1])

p_a(opacities; withlegend=false) = @pgf Axis(
    {
        "grid style" = {"color" = color_grid},
        "label style" = {"color" = color_text},
        "tick label style" = {"color" = color_text},
        "axis line style" = {"color" = color_axis},
        "legend style" = {"fill" = "none", "text" = color_text, "draw" = color_axis},
        label_style = raw"font=\LARGE",
        xmajorgrids,
        ymajorgrids,
        xlabel = L"$\Delta$ Semi-Major Axis (km)",
        ylabel = "RMS Position error (km)",
        legend_pos = "south east",
        ymode = "log"
    },
    PlotInc({lineopts..., opacity = opacities[1], color = color_CW}, Coordinates(Δa_range / 1000, pos_err_a_cw)),
    PlotInc({lineopts..., opacity = opacities[2], color = color_YA, style = "dashed"}, Coordinates(Δa_range / 1000, pos_err_a_ya)),
    PlotInc({lineopts..., opacity = opacities[3], color = color_KGD, style = "dashed"}, Coordinates(Δa_range / 1000, pos_err_a_kgd)),
    # PlotInc({lineopts..., opacity = opacities[4], color = color_LIN, style = "dotted"}, Coordinates(Δa_range / 1000, pos_err_a_lin)),
    PlotInc({lineopts..., opacity = opacities[5], color = color_KS}, Coordinates(Δa_range / 1000, pos_err_a_ks)),
    withlegend ? Legend(["CW", "YA", "KGD", "KS"]) : Legend()
)
p_a_cw = p_a([1, 0, 0, 0, 0])
p_a_ya = p_a([1, 1, 0, 0, 0])
p_a_kgd = p_a([1, 1, 1, 0, 0])
p_a_lin = p_a([1, 1, 1, 1, 0])
p_a_all = p_a([1, 1, 1, 1, 1])
p_a_all_legend = p_a([1, 1, 1, 1, 1]; withlegend=true)

pgfsave(joinpath("figs", "pdf", "rms_M_cw" * color_mode * ".pdf"), p_M_cw, dpi=300)
pgfsave(joinpath("figs", "pdf", "rms_M_ya" * color_mode * ".pdf"), p_M_ya, dpi=300)
pgfsave(joinpath("figs", "pdf", "rms_M_kgd" * color_mode * ".pdf"), p_M_kgd, dpi=300)
pgfsave(joinpath("figs", "pdf", "rms_M_lin" * color_mode * ".pdf"), p_M_lin, dpi=300)
pgfsave(joinpath("figs", "pdf", "rms_M_all" * color_mode * ".pdf"), p_M_all, dpi=300)
pgfsave(joinpath("figs", "pdf", "rms_i_cw" * color_mode * ".pdf"), p_i_cw, dpi=300)
pgfsave(joinpath("figs", "pdf", "rms_i_ya" * color_mode * ".pdf"), p_i_ya, dpi=300)
pgfsave(joinpath("figs", "pdf", "rms_i_kgd" * color_mode * ".pdf"), p_i_kgd, dpi=300)
pgfsave(joinpath("figs", "pdf", "rms_i_lin" * color_mode * ".pdf"), p_i_lin, dpi=300)
pgfsave(joinpath("figs", "pdf", "rms_i_all" * color_mode * ".pdf"), p_i_all, dpi=300)
pgfsave(joinpath("figs", "pdf", "rms_e_cw" * color_mode * ".pdf"), p_e_cw, dpi=300)
pgfsave(joinpath("figs", "pdf", "rms_e_ya" * color_mode * ".pdf"), p_e_ya, dpi=300)
pgfsave(joinpath("figs", "pdf", "rms_e_kgd" * color_mode * ".pdf"), p_e_kgd, dpi=300)
pgfsave(joinpath("figs", "pdf", "rms_e_lin" * color_mode * ".pdf"), p_e_lin, dpi=300)
pgfsave(joinpath("figs", "pdf", "rms_e_all" * color_mode * ".pdf"), p_e_all, dpi=300)
pgfsave(joinpath("figs", "pdf", "rms_a_cw" * color_mode * ".pdf"), p_a_cw, dpi=300)
pgfsave(joinpath("figs", "pdf", "rms_a_ya" * color_mode * ".pdf"), p_a_ya, dpi=300)
pgfsave(joinpath("figs", "pdf", "rms_a_kgd" * color_mode * ".pdf"), p_a_kgd, dpi=300)
pgfsave(joinpath("figs", "pdf", "rms_a_lin" * color_mode * ".pdf"), p_a_lin, dpi=300)
pgfsave(joinpath("figs", "pdf", "rms_a_all" * color_mode * ".pdf"), p_a_all, dpi=300)
pgfsave(joinpath("figs", "pdf", "rms_a_all_legend" * color_mode * ".pdf"), p_a_all_legend, dpi=300)
