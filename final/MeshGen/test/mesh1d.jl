#!/usr/bin/julia -f

import MeshGen: Vec, TileSet, Meshes, meshgen

immutable CircleModel{T} <: Meshes.Abstract1D{Vec{2,T}}
    scale::T
end

function Meshes.get_section{T}(model::CircleModel{T}, p::Vec{2,T})
    p.r[2] > 0 ? 1 : 0
end

function Meshes.get_step_size{T}(model::CircleModel{T},
                                 p::Vec{2,T}, section::Int)
    (p.r[2]^2 + T(0.2)) * model.scale * (1 + 1.5 * section)
end

function Meshes.get_init{T}(model::CircleModel{T})
    Vec{2,T}(0, -1), 0, (Vec{2,T}(1, 0),)
end

function Meshes.get_next_point{T}(model::CircleModel{T}, p::Vec{2,T},
                                  step::Vec{2,T}, section::Int, clip::Bool)
    next_p = step + p
    if next_p.r[2] * p.r[2] < 0
        section = section == 0 ? 1 : 0
        clip && (next_p = Vec{2,T}(next_p.r[1], 0))
    end
    next_p / abs(next_p), section
end

function Meshes.check_crossing{T}(model::CircleModel{T}, p1::Vec{2,T},
                                  p2::Vec{2,T}, tileset::TileSet{1,Vec{2,T}})
    pset = tileset.pts
    alg1 = atan2(p1.r[1], p1.r[2])
    alg2 = atan2(p2.r[1], p2.r[2])
    for (fidx, ftile) in tileset.tiles
        p = pset[ftile[1]]
        alg = atan2(p.r[1], p.r[2])
        diff1 = alg - alg1
        diff2 = alg - alg2
        while diff1 <= -π
            diff1 += π
        end
        while diff1 >= π
            diff1 -= π
        end
        while diff2 <= -π
            diff2 += π
        end
        while diff2 >= π
            diff2 -= π
        end
        diff1 * diff2 <= 0 && return Nullable((fidx, p))
    end
    return Nullable{Tuple{Int,Vec{2,T}}}()
end

ps = meshgen(CircleModel(0.1))

using PyPlot
pset = ps.pts
for (idx, tile) in ps.tiles
    r1, r2 = pset[tile]
    plot([r1.r[1], r2.r[1]], [r1.r[2], r2.r[2]], "g.-")
end
xlim([-1, 1])
ylim([-1, 1])
gca()[:set_aspect]("equal")
grid()
savefig("img/mesh1d.png", bbox_inches="tight", dpi=1000)
# show()
