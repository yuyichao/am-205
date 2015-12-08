#!/usr/bin/julia -f

module TestMeshes

import MeshGen: Vec, TileSet, Meshes, mesh

immutable CircleModel{T} <: Meshes.Abstract1D{Vec{2,T}}
    scale::T
end

function Meshes.get_section{T}(model::CircleModel{T}, p::Vec{2,T})
    p.r[2] > 0 ? 1 : 0
end

function Meshes.get_step_size{T}(model::CircleModel{T},
                                 p::Vec{2,T}, section::Int)
    (p.r[2]^2 + T(0.1)) * model.scale * (1 + 2 * section)
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

ps = mesh(CircleModel(0.1))

# using PyPlot
# pset = ps.pts
# for (idx, tile) in ps.tiles
#     r1, r2 = pset[tile]
#     plot([r1.r[1], r2.r[1]], [r1.r[2], r2.r[2]], "bo-")
# end
# # plot(sin(ps), cos(ps), "o")
# grid()
# show()

immutable CircleModel2D{T} <: Meshes.Abstract2D{Vec{2,T}}
end

function Meshes.get_section{T}(::CircleModel2D{T}, p::Vec{2,T})
    abs(p) >= 1 ? -1 : 0
end
function Meshes.get_step_size{T}(::CircleModel2D{T}, p::Vec{2,T}, section::Int)
    T(0.05) / T(abs(p) + 0.3)
end

function Meshes.get_init{T}(::CircleModel2D{T})
    Vec{2,T}(0, 0), 0, (Vec{2,T}(1, 0), Vec{2,T}(0.5, √(3) / 2))
end

function Meshes.get_next_point{T}(::CircleModel2D{T}, p::Vec{2,T},
                                  step::Vec{2,T}, section::Int, clip::Bool)
    next_p = p + step
    section = abs(next_p) >= 1 ? -1 : 0
    if section == -1 && clip
        next_p /= abs(next_p)
    end
    next_p, section
end

function same_side_of(lp1::Vec{2}, lp2::Vec{2}, p1::Vec{2}, p2::Vec{2})
    # http://stackoverflow.com/a/7069702/904262
    lv = lp2 - lp1
    p11 = p1 - lp1
    p21 = p2 - lp1
    # which side is p1 on calculate lv × p11
    side1 = lv.r[1] * p11.r[2] - lv.r[2] * p11.r[1]
    # which side is p2 on calculate lv × p21
    side2 = lv.r[1] * p21.r[2] - lv.r[2] * p21.r[1]
    side1 * side2 > 0
end

function line_crosses(l1p1, l1p2, l2p1, l2p2)
    # http://stackoverflow.com/a/7069702/904262
    same_side_of(l1p1, l1p2, l2p1, l2p2) && return false
    same_side_of(l2p1, l2p2, l1p1, l1p2) && return false
    return true
end

function point_in_angle(common, p1, p2, ptest)
    same_side_of(common, ptest, p1, p2) && return false
    v1 = p1 - common
    v2 = p2 - common
    vtest = ptest - common
    v1 /= abs(v1)
    v2 /= abs(v2)
    vtest /= abs(vtest)
    v1 * vtest <= 0 && v2 * vtest <= 0 && return false
    acos(v1 * vtest) + acos(v2 * vtest) <= π
end

function Meshes.check_crossing{T}(::CircleModel2D{T}, ps::NTuple{2,Vec{2,T}},
                                  p2::Vec{2,T}, tileset::TileSet{2,Vec{2,T}})
    crossing = Int[]
    pset = tileset.pts
    pidx1 = pset[ps[1]]
    pidx2 = pset[ps[2]]
    for (fidx, ftile) in tileset.tiles
        ftidx1, ftidx2 = ftile
        if ((ftidx1 == pidx1 && ftidx2 == pidx2) ||
            (ftidx1 == pidx2 && ftidx2 == pidx1))
            continue
        end
        ftp1, ftp2 = pset[ftile]
        if ((ftp1 == ps[1] && ftp2 == p2) || (ftp1 == p2 && ftp2 == ps[1]))
            continue
        end
        if ((ftp1 == ps[2] && ftp2 == p2) || (ftp1 == p2 && ftp2 == ps[2]))
            continue
        end
        if pidx1 == ftidx1
            point_in_angle(ps[1], p2, ps[2], ftp2) && push!(crossing, fidx)
            continue
        elseif pidx1 == ftidx2
            point_in_angle(ps[1], p2, ps[2], ftp1) && push!(crossing, fidx)
            continue
        elseif pidx2 == ftidx1
            point_in_angle(ps[2], p2, ps[1], ftp2) && push!(crossing, fidx)
            continue
        elseif pidx2 == ftidx2
            point_in_angle(ps[2], p2, ps[1], ftp1) && push!(crossing, fidx)
            continue
        elseif ftp1 == p2
            point_in_angle(p2, ps[2], ps[1], ftp2) && push!(crossing, fidx)
            continue
        elseif ftp2 == p2
            point_in_angle(p2, ps[2], ps[1], ftp1) && push!(crossing, fidx)
            continue
        end
        if (line_crosses(ftp1, ftp2, ps[1], p2) ||
            line_crosses(ftp1, ftp2, ps[2], p2))
            push!(crossing, fidx)
        end
    end
    crossing
end

ps2d, frontier = mesh(CircleModel2D{Float64}())

using PyPlot
pset2 = ps2d.pts
figure()
for (idx, tile) in ps2d.tiles
    r1, r2, r3 = pset2[tile]
    plot([r1.r[1], r2.r[1], r3.r[1], r1.r[1]],
         [r1.r[2], r2.r[2], r3.r[2], r1.r[2]], "o-")
end
grid()
figure()
for (idx, tile) in frontier.tiles
    r1, r2 = pset2[tile]
    plot([r1.r[1], r2.r[1]], [r1.r[2], r2.r[2]], "bo-")
    # if idx in ignored
    #     plot([r1.r[1], r2.r[1]], [r1.r[2], r2.r[2]], "go-")
    # else
    # end
end
grid()
show()

end
