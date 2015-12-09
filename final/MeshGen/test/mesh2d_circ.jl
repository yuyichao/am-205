#!/usr/bin/julia -f

import MeshGen: Vec, TileSet, Meshes, meshgen

immutable CircleWithHoleModel2D{T} <: Meshes.Abstract2D{Vec{2,T}}
    hole_c::Vec{2,T}
    hole_r::T
end

function Meshes.get_section{T}(m::CircleWithHoleModel2D{T}, p::Vec{2,T})
    if abs(p) >= 1
        return -1
    end
    if abs(p - m.hole_c) <= m.hole_r
        return -1
    end
    return 0
end
function Meshes.get_step_size{T}(m::CircleWithHoleModel2D{T}, p::Vec{2,T},
                                 section::Int)
    density = abs(p) + sqrt(m.hole_r / abs(p - m.hole_c)^3) + 0.3
    T(0.1) / T(density)
end

function Meshes.get_init{T}(m::CircleWithHoleModel2D{T})
    init_p = m.hole_c / abs(m.hole_c) * (abs(m.hole_c) + m.hole_r + 1) / 2
    init_p, 0, (Vec{2,T}(1, 0), Vec{2,T}(0.5, √(3) / 2))
end

function Meshes.get_next_point{T}(m::CircleWithHoleModel2D{T}, p::Vec{2,T},
                                  step::Vec{2,T}, section::Int, clip::Bool)
    next_p = p + step
    if abs(next_p) >= 1
        section = -1
        clip && (next_p /= abs(next_p))
    elseif abs(next_p - m.hole_c) <= m.hole_r
        section = -1
        if clip
            hole_p = next_p - m.hole_c
            next_p = m.hole_c + hole_p / abs(hole_p) * m.hole_r
        end
    else
        section = 0
    end
    next_p, section
end

function same_side_of(lp1::Vec{2}, lp2::Vec{2}, p1::Vec{2}, p2::Vec{2})
    # http://stackoverflow.com/a/7069702/904262
    lv = lp2 - lp1
    p11 = p1 - lp1
    p21 = p2 - lp1
    # which side is p1 on
    side1 = lv × p11
    # which side is p2 on
    side2 = lv × p21
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

function Meshes.check_crossing{T}(m::CircleWithHoleModel2D{T},
                                  ps::NTuple{2,Vec{2,T}},
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


@time ps2d, frontier =
    meshgen(CircleWithHoleModel2D{Float64}(Vec{2,Float64}(0.0, 0.6), 0.2))

using PyPlot
pset2 = ps2d.pts
# figure()
# for (idx, tile) in ps2d.tiles
#     r1, r2, r3 = pset2[tile]
#     plot([r1.r[1], r2.r[1], r3.r[1], r1.r[1]],
#          [r1.r[2], r2.r[2], r3.r[2], r1.r[2]], "-")
# end
# grid()
# xlim([-1, 1])
# ylim([-1, 1])
# gca()[:set_aspect]("equal")
# savefig("mesh2d_circ.png", bbox_inches="tight")

figure()
for (idx, tile) in frontier.tiles
    r1, r2 = pset2[tile]
    # ewidth = sqrt(abs(r1 - r2)) * 2
    ewidth = 1
    plot([r1.r[1], r2.r[1]], [r1.r[2], r2.r[2]], linewidth=ewidth, "g-")
end
grid()
xlim([-0.25, -0.03])
ylim([0.05, 0.17])
gca()[:set_aspect]("equal")
savefig("mesh2d_circ_opt.png", bbox_inches="tight", dpi=1000)
show()
