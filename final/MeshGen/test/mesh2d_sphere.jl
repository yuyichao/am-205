#!/usr/bin/julia -f

import MeshGen: Vec, TileSet, Meshes, meshgen

immutable SphereModel2D{T} <: Meshes.Abstract2D{Vec{3,T}}
end

Meshes.get_section{T}(m::SphereModel2D{T}, p::Vec{3,T}) = 0
Meshes.get_step_size{T}(m::SphereModel2D{T}, p::Vec{3,T}, section::Int) =
    T(0.15)

function Meshes.get_init{T}(m::SphereModel2D{T})
    Vec{3,T}(0, 0, -1), 0, (Vec{3,T}(1, 0, 0), Vec{3,T}(0, 1, 0))
end

function Meshes.get_next_point{T}(m::SphereModel2D{T}, p::Vec{3,T},
                                  step::Vec{3,T}, section::Int, clip::Bool)
    next_p = p + step
    next_p / abs(next_p), section
end

function same_side_of(lp1::Vec{3}, lp2::Vec{3}, p1::Vec{3}, p2::Vec{3})
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

function Meshes.check_crossing{T}(m::SphereModel2D{T},
                                  ps::NTuple{2,Vec{3,T}},
                                  p2::Vec{3,T}, tileset::TileSet{2,Vec{3,T}})
    crossing = Int[]
    pset = tileset.pts
    pidx1 = pset[ps[1]]
    pidx2 = pset[ps[2]]
    for (fidx, ftile) in tileset.tiles
        ftidx1, ftidx2 = ftile
        # Ignore the edges we are working on
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
        # handle the case where we have a shared point
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


@time ps2d, frontier = meshgen(SphereModel2D{Float64}())

using PyPlot

pset2 = ps2d.pts
figure()
for (idx, tile) in ps2d.tiles
    r1, r2, r3 = pset2[tile]
    plot3D([r1.r[1], r2.r[1], r3.r[1], r1.r[1]],
           [r1.r[2], r2.r[2], r3.r[2], r1.r[2]],
           [r1.r[3], r2.r[3], r3.r[3], r1.r[3]], "-")
end
grid()
figure()
for (idx, tile) in frontier.tiles
    r1, r2 = pset2[tile]
    plot3D([r1.r[1], r2.r[1]], [r1.r[2], r2.r[2]], [r1.r[3], r2.r[3]], "b-")
end
grid()

figure()
for (idx, tile) in ps2d.tiles
    r1, r2, r3 = pset2[tile]
    (r1.r[3] >= 0 && r2.r[3] >= 0 && r3.r[3] >= 0) || continue
    plot([r1.r[1], r2.r[1], r3.r[1], r1.r[1]],
         [r1.r[2], r2.r[2], r3.r[2], r1.r[2]], "-")
end
grid()
figure()
for (idx, tile) in ps2d.tiles
    r1, r2, r3 = pset2[tile]
    (r1.r[3] < 0 || r2.r[3] < 0 && r3.r[3] < 0) || continue
    plot([r1.r[1], r2.r[1], r3.r[1], r1.r[1]],
         [r1.r[2], r2.r[2], r3.r[2], r1.r[2]], "-")
end
grid()

figure()
for (idx, tile) in frontier.tiles
    r1, r2 = pset2[tile]
    (r1.r[3] >= 0 && r2.r[3] >= 0) || continue
    plot([r1.r[1], r2.r[1]], [r1.r[2], r2.r[2]], "b-")
end
grid()
figure()
for (idx, tile) in frontier.tiles
    r1, r2 = pset2[tile]
    (r1.r[3] < 0 || r2.r[3] < 0) || continue
    plot([r1.r[1], r2.r[1]], [r1.r[2], r2.r[2]], "b-")
end
grid()

show()
