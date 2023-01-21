using Pkg;
Pkg.activate(joinpath(@__DIR__, ".."));
using JLD2
using PGFPlotsX
using LaTeXStrings
using Colors
using LinearAlgebra

SAVEAS_PNG = true

color_list = distinguishable_colors(6, [RGB(1, 1, 1), RGB(0, 0, 0)], dropseed=true)
const def_linewidth = "ultra thick"
const color_fixed = color_list[4]
const color_newton = color_list[2]
# lineopts = @pgf {no_marks, "very thick", "error bars/y dir=both", "error bars/y explicit"}
lineopts = @pgf {no_marks, "ultra thick"}

data = load("data/cart_to_ks_transform_smoothness.jld2");


times = data["times"]
x_traj_ks_fixed = data["x_traj_ks_fixed"]
x_traj_ks_newton = data["x_traj_ks_newton"]

x_traj_ks_fixed = hcat(x_traj_ks_fixed...)
x_traj_ks_newton = hcat(x_traj_ks_newton...)

# put things into a way pgf likes
combined_times = [times; NaN; times; NaN; times; NaN; times]
combined_fixed = [x_traj_ks_fixed[1, :]; NaN; x_traj_ks_fixed[2, :]; NaN; x_traj_ks_fixed[3, :]; NaN; x_traj_ks_fixed[4, :]]
combined_newton = [x_traj_ks_newton[1, :]; NaN; x_traj_ks_newton[2, :]; NaN; x_traj_ks_newton[3, :]; NaN; x_traj_ks_newton[4, :]]

p = @pgf Axis(
    {
        xmajorgrids,
        ymajorgrids,
        xlabel = "Time (nondimensional)",
        ylabel = "KS Position",
        legend_pos = "south west",
    },
    PlotInc({no_marks, "ultra thick", color = color_fixed}, Coordinates(combined_times, combined_fixed)),
    PlotInc({no_marks, "ultra thick", color = color_newton, style = "dotted"}, Coordinates(combined_times, combined_newton)),
    Legend(["Common", "Nearest"])
)

if SAVEAS_PNG
    pgfsave(joinpath("figs", "png", "cart_to_ks_transform_smoothness.png"), p, dpi=300)
else
    pgfsave(joinpath("figs", "cart_to_ks_transform_smoothness.tikz"), p, include_preamble=false)
end