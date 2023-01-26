using MeshCat
using GeometryBasics
using Meshes
using CoordinateTransformations, Rotations
using Colors
using JLD2
using LinearAlgebra

model_scale = 0.005

chaser_material = MeshPhongMaterial(wireframe=false, wireframeLinewidth=0.001, side=2, color=RGBA(0.2, 0.9, 0.9, 0.8))
target_material = MeshPhongMaterial(wireframe=false, wireframeLinewidth=0.001, side=2, color=RGBA(0.9, 0.9, 0.2, 0.8))

vis = Visualizer()
render(vis)

setprop!(vis["/Background"], "visible", false)
setprop!(vis["/Grid"], "visible", false)

satobj = MeshFileGeometry("anim/sat_model.stl")

rotation = RotZ(-pi / 2) * RotY(-pi / 4)
scale_matrix = I(3) * model_scale

setobject!(vis["chaser"], satobj, chaser_material)
settransform!(vis["chaser"], compose(Translation(0, 0, 0), LinearMap(rotation)))
setprop!(vis["chaser"], "scale", model_scale * ones(3))

setobject!(vis["target"], satobj, target_material)
settransform!(vis["target"], compose(Translation(0, 0, 0), LinearMap(rotation)))
setprop!(vis["target"], "scale", model_scale * ones(3))


data = load("data/low_thrust_rendezvous.jld2");
scales = data["scales"]

d_scale = scales.d_scale
t_scale = scales.t_scale

relative_states_rtn = data["relative_states_rtn"]
times = data["times"]
Ndata = length(times)

anim = Animation(120)

for i = 1:Ndata
    R = 0.1 * (d_scale / 1000) * relative_states_rtn[i][1]
    T = 0.1 * (d_scale / 1000) * relative_states_rtn[i][2]
    N = 0.1 * (d_scale / 1000) * relative_states_rtn[i][3]

    atframe(anim, i) do
        settransform!(vis["chaser"], compose(Translation(T, N, R), LinearMap(rotation)))
        setprop!(vis["chaser"], "scale", model_scale * ones(3))
    end
end

setanimation!(vis, anim)


delete!(vis)

####
# To save Animation, use  Animations > Recording > record in the viewer
# To convert animation frames to mp4, 
# ffmpeg -r 60 -i %07d.png -vcodec libx264 -pix_fmt yuv420p -preset slow -crf 18 ../rendezvous_animation.mp4