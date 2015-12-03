#!/usr/bin/julia -f

using MeshGen

immutable CircleModel{T} <: MeshGen.Abstract1D{T}
    scale::T
end

MeshGen.get_section{T}(m::CircleModel{T}, p) =
    mod(floor(Int, p / π + T(0.5)), 2)
MeshGen.get_step_size{T}(model::CircleModel{T}, p, cur_section) =
    (sin(p)^2 + T(0.1)) * model.scale * (1 + 4 * cur_section)
MeshGen.get_edge{T}(model::CircleModel{T}, cur_section,
                    next_section, p, next_p) =
    π * (floor(next_p / π + T(0.5)) - T(0.5))
MeshGen.check_crossing(model, p, frontier) =
    p < frontier[2] || p > frontier[1] + 2π

ps = MeshGen.mesh_1d(CircleModel(0.1))

# using PyPlot
# plot(sin(ps), cos(ps), "o")
# grid()
# show()
