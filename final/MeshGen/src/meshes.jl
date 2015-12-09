#!/usr/bin/julia -f

module Meshes

import ..TileSets: PointSet, TileSet

abstract AbstractModel{N,V}

"""
    get_section{N,V}(model::AbstractModel{N,V}, point::V) -> Int

When the `model` has multiple sections, returns which section the `point`
belongs to. Return `-1` for out of bound points.
"""
function get_section end

"""
    get_step_size{N,V}(model::AbstractModel{N,V},
                       point::V, cur_section::Int)

Returns the step size for the current `point`, When the `point` is on the
boundary between two sections with different step size, `cur_section` determines
which section is used. Return type is the same with the element type of `V`.
"""
function get_step_size end

# """
#     get_bases{N,V}(model::AbstractModel{N,V}, point::V) -> NTuple{N,V}

# Returns a set of base vectors in the `model` geometry at `point`.
# """
# function get_bases end

"""
    get_init{N,V}(model::AbstractModel{N,V}) -> Tuple{V,Int,NTuple{N,V}}

Returns parameters to initialize the meshing process, including the coordinate
of the first point, the section it belongs to and the initial directions to
move.
"""
function get_init end

"""
    get_next_point{N,V}(model::AbstractModel{N,V}, point::V,
                        step::V, section::Int, clip::Bool) -> Tuple{V,Int}

Returns the coordinate and the section of the next point that is obtain by
moving from `point` by `step` (vector). If `clip` is `true` and the line
crosses a boundary, the returned coordinate should be on the boundary.
If there's multiple boundaries, return the first crossing.
"""
function get_next_point end

"""
    check_crossing{N,V}(model::AbstractModel{N,V}, point1::V, point2::V,
                        tileset::TileSet{N,V}) -> Nullable{Tuple{Int,V}}

Check if the line from `point1` to `point2` intersect with any tile in `tileset`
Return a null `Nullable` if not. Otherwise, return the index of the intersecting
tile and the coordinate of the intersection
"""
function check_crossing end

"""
    check_crossing{N,V}(model::AbstractModel{N,V}, orig_points::NTuple{N,V},
                        point2::V, tileset::TileSet{N,V})
        -> Vector{Int}

Check if the lines from `orig_point` to `point2` intersect with any tile in
`tileset`. Return an empty `Vector` if not. Otherwise, return the indecis of
the intersecting tiles.
"""
function check_crossing end

typealias Abstract1D{V} AbstractModel{1,V}

function init_mesh{V}(model::Abstract1D{V}, mesh, frontier, dirs, sections)
    # Get initial point
    p, section, (dir,) = get_init(model)
    step = get_step_size(model, p, section)
    # try move in both directions
    next_p1, next_sec1 = get_next_point(model, p, step * dir, section, true)
    next_p2, next_sec2 = get_next_point(model, p, -step * dir, section, true)
    if next_sec1 != -1
        push!(mesh, (p, next_p1))
        fidx = push!(frontier, (next_p1,))
        dirs[fidx] = dir
        sections[fidx] = next_sec1
    end
    if next_sec2 != -1
        push!(mesh, (p, next_p2))
        fidx = push!(frontier, (next_p2,))
        dirs[fidx] = -dir
        sections[fidx] = next_sec2
    end
    nothing
end

function meshgen{V}(model::Abstract1D{V})
    pset = PointSet{V}()
    mesh = TileSet{2}(pset)
    frontier = TileSet{1}(pset)
    dirs = Dict{Int,V}()
    sections = Dict{Int,Int}()
    init_mesh(model, mesh, frontier, dirs, sections)

    while !isempty(frontier.tiles)
        fidx, ftile = first(frontier.tiles)
        dir = dirs[fidx]
        section = sections[fidx]
        delete!(frontier, fidx)
        delete!(dirs, fidx)
        delete!(sections, fidx)
        section == -1 && continue

        p = pset[ftile[1]]
        step = get_step_size(model, p, section)
        @assert step > 0

        # First go twice as far to figure out what might be closed to the
        # next point.
        next_p, next_sec = get_next_point(model, p, step * 2 * dir,
                                          section, false)
        # Now check if this line intersect with the current frontier
        iscross = check_crossing(model, p, next_p, frontier)
        if isnull(iscross)
            # Not cross
            # Now check if there's a boundary nearby
            if next_sec == section
                # No there isn't. Now we can simply add the new point at step
                next_p, next_sec = get_next_point(model, p, step * dir,
                                                  section, true)
                @assert next_sec == section
                diff_vec = next_p - p
                push!(mesh, (p, next_p))
                next_fidx = push!(frontier, (next_p,))
                diff_len = abs(diff_vec)
                @assert diff_len > 0
                dirs[next_fidx] = diff_vec / diff_len
                sections[next_fidx] = section
            else
                # OK we hit a boundary, figure out how far it is
                next_p, next_sec2 = get_next_point(model, p, step * 2 * dir,
                                                   section, true)
                @assert next_sec2 == next_sec
                diff_vec = next_p - p
                step_len = abs(diff_vec)
                step_len <= step / 10 && continue
                if step_len < step * √(2)
                    # if the boundary is closed enough, simply use that point
                    push!(mesh, (p, next_p))
                    next_fidx = push!(frontier, (next_p,))
                    dirs[next_fidx] = diff_vec / step_len
                    sections[next_fidx] = next_sec
                else
                    # otherwise, try the middle point.
                    next_p, next_sec = get_next_point(model, p,
                                                      step_len * dir / 2,
                                                      section, true)
                    @assert next_sec == section
                    diff_vec = next_p - p
                    step_len = abs(diff_vec)
                    @assert step_len > 0
                    push!(mesh, (p, next_p))
                    next_fidx = push!(frontier, (next_p,))
                    dirs[next_fidx] = diff_vec / step_len
                    sections[next_fidx] = section
                end
            end
        else
            # cross
            cross_tidx, cross_p = get(iscross)
            # Now check if we've also crossed a boundary
            if next_sec == section
                # No there isn't. We only need to coordinate the two frontiers
                diff_vec = cross_p - p
                step_len = abs(diff_vec)
                @assert step_len > 0
                if step_len < step * √(2)
                    # if it's close enough, simply connect the two
                    push!(mesh, (p, cross_p))
                    delete!(frontier, cross_tidx)
                    delete!(dirs, cross_tidx)
                    delete!(sections, cross_tidx)
                else
                    # otherwise, split into two and connect them
                    next_p, next_sec = get_next_point(model, p,
                                                      step_len / 2 * dir,
                                                      section, true)
                    @assert next_sec == section
                    push!(mesh, (p, next_p))
                    push!(mesh, (next_p, cross_p))
                    delete!(frontier, cross_tidx)
                    delete!(dirs, cross_tidx)
                    delete!(sections, cross_tidx)
                end
            else
                # So we've also crossed a boundary, let's see where we crossed
                # it
                next_p, next_sec2 = get_next_point(model, p, step * 2 * dir,
                                                   section, true)
                @assert next_sec == next_sec2
                next_len = abs(next_p - p)
                cross_len = abs(cross_p - p)
                @assert next_len > 0
                @assert cross_len > 0
                if cross_len <= next_len * √(2)
                    # if they are closed enough or the crossing with the
                    # frontier is closer, just connect the next point
                    push!(mesh, (p, cross_p))
                    delete!(frontier, cross_tidx)
                    delete!(dirs, cross_tidx)
                    delete!(sections, cross_tidx)
                else
                    # We are closer to the boundary, connect the boundary and
                    # then the frontier
                    push!(mesh, (p, next_p))
                    push!(mesh, (next_p, cross_p))
                    delete!(frontier, cross_tidx)
                    delete!(dirs, cross_tidx)
                    delete!(sections, cross_tidx)
                end
            end
        end
    end

    mesh
end

typealias Abstract2D{T} AbstractModel{2,T}

immutable WorkSet2D{V}
    edges::TileSet{2,V} # All edges
    fr::Set{Int} # Current frontier
    dirs::Dict{Int,V} # Directions
    psect::Dict{Int,Set{Int}} # Sections of points
    esect::Dict{Int,Int} # Sections of points
    p2e::Dict{Int,Set{Int}} # Map from point to edges
    WorkSet2D(pset) = new(TileSet{2}(pset), Set{Int}(), Dict{Int,V}(),
                          Dict{Int,Set{Int}}(), Dict{NTuple{2,Int},Int}(),
                          Dict{Int,Set{Int}}())
end

function get_psects!{V}(ws::WorkSet2D{V}, pidx::Int)
    if pidx in keys(ws.psect)
        ws.psect[pidx]
    else
        ws.psect[pidx] = Set{Int}()
    end
end

function add_edge_section!{V}(ws::WorkSet2D{V}, eidx::Int, section::Int)
    ws.esect[eidx] = section
end

function add_section!{V}(ws::WorkSet2D{V}, pidx::Int, section::Int)
    sects = get_psects!(ws, pidx)
    push!(sects, section)
    nothing
end

add_section!{V}(ws::WorkSet2D{V}, p::V, section::Int) =
    add_section!(ws, ws.edges.pts[p], section)

function get_section{V}(ws::WorkSet2D{V}, pidx1::Int, pidx2::Int)
    eidx = ws.edges[(pidx1, pidx2)]
    sects1 = get_psects!(ws, pidx1)
    sects2 = get_psects!(ws, pidx2)
    if eidx in keys(ws.esect)
        section = ws.esect[eidx]
        sects1 = delete!(copy(sects1), section)
        sects2 = delete!(copy(sects2), section)
        sec1 = isempty(sects1) ? section : first(sects1)
        sec2 = isempty(sects2) ? section : first(sects2)
        return sec1, sec2, section
    end
    if isempty(sects1)
        isempty(sects2) && return (-1, -1, -1)
        sec = first(sects2)
        return (-1, sec, sec)
    elseif isempty(sects2)
        sec = first(sects1)
        return (sec, -1, sec)
    end
    for sec in sects1
        sec in sects2 && return (sec, sec, sec)
    end
    sec1 = first(sects1)
    return (sec1, first(sects2), sec1)
end
get_section{V}(ws::WorkSet2D{V}, p1::V, p2::V) =
    get_section(ws, ws.edges.pts[p1], ws.edges.pts[p2])

function Base.push!{V}(ws::WorkSet2D{V}, pidxs::NTuple{2,Int}, dir::V)
    eidx = push!(ws.edges, pidxs)
    push!(ws.fr, eidx)
    ws.dirs[eidx] = dir
    for pidx in pidxs
        push!(get!(ws.p2e, pidx, Set{Int}()), eidx)
    end
    eidx
end
Base.push!{V}(ws::WorkSet2D{V}, ps::NTuple{2,V}, dir::V) =
    push!(ws, ws.edges.pts[ps], dir)

Base.isempty(ws::WorkSet2D) = isempty(ws.fr)

Base.start(ws::WorkSet2D) = start(ws.fr)
function Base.next(ws::WorkSet2D, idx)
    eidx, idx = next(ws.fr, idx)
    pidxs = ws.edges[eidx]
    dir = ws.dirs[eidx]
    (eidx, ws.edges.pts[pidxs], dir), idx
end
Base.done(ws::WorkSet2D, idx) = done(ws.fr, idx)

function Base.delete!(ws::WorkSet2D, eidx::Int)
    delete!(ws.fr, eidx)
    delete!(ws.dirs, eidx)
end
function Base.getindex(ws::WorkSet2D, eidx::Int)
    pidxs = ws.edges[eidx]
    dir = ws.dirs[eidx]
    ws.edges.pts[pidxs], dir
end

function init_mesh{V}(model::Abstract2D{V}, mesh, ws)
    # Get initial point
    p, section, (dir1, dir2) = get_init(model)
    step = get_step_size(model, p, section)
    next_p1, next_sec1 = get_next_point(model, p, step * dir1, section, true)
    next_p2, next_sec2 = get_next_point(model, p, step * dir2, section, true)
    push!(mesh, (p, next_p1, next_p2))
    dir3 = dir1 + dir2
    dir3 /= -abs(dir3)
    dir1, dir2 = calc_directions(next_p1, next_p2, p)
    eidx1 = push!(ws, (next_p1, p), dir1)
    eidx2 = push!(ws, (next_p2, p), dir2)
    eidx12 = push!(ws, (next_p1, next_p2), dir3)
    add_section!(ws, next_p1, next_sec1)
    add_section!(ws, next_p2, next_sec2)
    add_section!(ws, p, section)
    add_section!(ws, next_p1, section)
    add_section!(ws, next_p2, section)
    add_edge_section!(ws, eidx1, section)
    add_edge_section!(ws, eidx2, section)
    add_edge_section!(ws, eidx12, section)
    nothing
end

function calc_directions(old_p1, old_p2, new_p)
    edge1 = new_p - old_p1
    edge2 = new_p - old_p2
    ledge1 = abs(edge1)
    ledge2 = abs(edge2)
    @assert ledge1 > 0 && ledge2 > 0
    edge1 /= ledge1
    edge2 /= ledge2

    proj = edge1 * edge2

    dir1 = edge2 - edge1 * proj
    dir2 = edge1 - edge2 * proj
    ldir1 = abs(dir1)
    ldir2 = abs(dir2)
    @assert ldir1 > 0 && ldir2 > 0
    dir1 / ldir1, dir2 / ldir2
end

function normalized_dot(v1, v2)
    v1 /= abs(v1)
    v2 /= abs(v2)
    v1 * v2
end

function orthorg_vec(v1, v2)
    n1 = v1 / abs(v1)
    v2 - (n1 * v2) * n1
end

function distance_to_segment(v1, v2)
    # p0 -> p1: v1
    # p0 -> p2: v2
    # p1 -> p2: v2 - v1 = -v12
    v12 = v1 - v2

    # check if angle 1 is larger than 90deg
    if v1 * v12 <= 0
        return abs(v1)
    end
    # check if angle 2 is larger than 90deg
    if v2 * v12 >= 0
        return abs(v2)
    end
    # use the distance
    abs(orthorg_vec(v12, v1))
end

function check_neighbers(model, ws, pset, mesh, p1, p2, dir, section)
    pidx1 = pset[p1]
    pidx2 = pset[p2]
    for eidx in union(ws.p2e[pidx1], ws.p2e[pidx2])
        etile = ws.edges.tiles[eidx]
        ep1, ep2 = pset[etile]
        ((ep1 == p1 && ep2 == p2) || (ep1 == p2 && ep2 == p1)) && continue
        if ep1 == p1
            p_com = p1
            p_old = p2
            p_new = ep2
        elseif ep2 == p1
            p_com = p1
            p_old = p2
            p_new = ep1
        elseif ep1 == p2
            p_com = p2
            p_old = p1
            p_new = ep2
        elseif ep2 == p2
            p_com = p2
            p_old = p1
            p_new = ep1
        else
            continue
        end
        (p_new - p_com) * dir <= 0 && continue
        normalized_dot(p_new - p_com, p_old - p_com) < 0.5 && continue
        isempty(check_crossing(model, (p_old, p_com), p_new,
                               ws.edges)) || continue
        push!(mesh, (p_old, p_new, p_com))
        # delete the other edge
        delete!(ws, eidx)
        try
            # Check if the thrid edge also exist
            eidx2 = ws.edges[(p_old, p_new)]
            delete!(ws, eidx2)
            return true
        end
        dir_pold, dir_pnew = calc_directions(p_com, p_old, p_new)
        eidx_new = push!(ws, (p_old, p_new), dir_pnew)
        add_section!(ws, p_new, section)
        add_edge_section!(ws, eidx_new, section)
        return true
    end
    return false
end

function check_free_2x(model, ws, pset, mesh, p1, p2, step, dir, section)
    pmid = (p1 + p2) / 2
    p_new, next_sec = get_next_point(model, pmid, step * dir * 2,
                                     section, false)
    # Give up if there's any crossing
    isempty(check_crossing(model, (p1, p2), p_new, ws.edges)) || return false
    if abs(p_new - pmid) >= step
        p_new, next_sec = get_next_point(model, pmid, (p_new - pmid) / 2,
                                         section, true)
    else
        p_new, next_sec = get_next_point(model, pmid,
                                         (p_new - pmid) * oftype(step, 0.9),
                                         section, true)
    end
    if abs(p_new - pmid) < 0.3 * step || next_sec != section
        return false
    end
    push!(mesh, (p1, p2, p_new))
    dir_p1, dir_p2 = calc_directions(p1, p2, p_new)
    eidx_new1 = push!(ws, (p_new, p1), dir_p1)
    eidx_new2 = push!(ws, (p_new, p2), dir_p2)

    add_section!(ws, p_new, next_sec)
    add_edge_section!(ws, eidx_new1, section)
    add_edge_section!(ws, eidx_new2, section)
    return true
end

# const do_abort = Ref(false)

function check_in_circle(model, ws, pset, mesh, p1, p2, step, dir, section)
    pmid = (p1 + p2) / 2
    dir = orthorg_vec(p1 - p2, dir)
    dir = dir / abs(dir)
    po = pmid + step * dir * oftype(step, 0.8)
    ro = abs(po - p1)
    remset = Set{Int}()
    # FIXME
    for (eidx, etile) in ws.edges.tiles
        ep1, ep2 = pset[etile]
        ((ep1 == p1 && ep2 == p2) || (ep1 == p2 && ep2 == p1)) && continue
        # Check if the edge is in the circle
        distance_to_segment(ep1 - po, ep2 - po) >= ro && continue
        ep1_valid = (ep1 - pmid) * dir > 0 && !(ep1 == p1 || ep1 == p2)
        ep2_valid = (ep2 - pmid) * dir > 0 && !(ep2 == p1 || ep2 == p2)
        cross1 = if ep1_valid
            check_crossing(model, (p1, p2), ep1, ws.edges)
        else
            Int[]
        end
        cross2 = if ep2_valid
            check_crossing(model, (p1, p2), ep2, ws.edges)
        else
            Int[]
        end
        ep1_valid &= isempty(cross1)
        ep2_valid &= isempty(cross2)
        if !(ep1_valid || ep2_valid)
            union!(remset, cross1)
            union!(remset, cross2)
            continue
        end
        p_new = ep1_valid ? ep1 : ep2

        push!(mesh, (p1, p2, p_new))
        # delete other edges
        eg1_exists = false
        eg2_exists = false
        try
            eidx1 = ws.edges[(p1, p_new)]
            delete!(ws, eidx1)
            eg1_exists = true
        end
        try
            eidx2 = ws.edges[(p2, p_new)]
            delete!(ws, eidx2)
            eg2_exists = true
        end
        dir_p1, dir_p2 = calc_directions(p1, p2, p_new)
        if !eg1_exists
            eidx_new1 = push!(ws, (p1, p_new), dir_p1)
            add_edge_section!(ws, eidx_new1, section)
        end
        if !eg2_exists
            eidx_new2 = push!(ws, (p2, p_new), dir_p2)
            add_edge_section!(ws, eidx_new2, section)
        end
        add_section!(ws, p_new, section)
        return true
    end
    isempty(remset) && return false
    threshold = oftype(step, Inf)
    remset2 = Set{Int}()
    while !isempty(remset)
        new_thresh = oftype(step, -Inf)
        for eidx in remset
            epidx1, epidx2 = ws.edges[eidx]
            ep1 = pset[epidx1]
            ep2 = pset[epidx2]
            ep1_valid = (ep1 - pmid) * dir > 0 && !(ep1 == p1 || ep1 == p2)
            ep2_valid = (ep2 - pmid) * dir > 0 && !(ep2 == p1 || ep2 == p2)
            dot1 = normalized_dot(ep1 - p1, ep1 - p2)
            dot2 = normalized_dot(ep2 - p1, ep2 - p2)
            ep1_valid &= dot1 < threshold
            ep2_valid &= dot2 < threshold

            cross1 = if ep1_valid
                check_crossing(model, (p1, p2), ep1, ws.edges)
            else
                Int[]
            end
            cross2 = if ep2_valid
                check_crossing(model, (p1, p2), ep2, ws.edges)
            else
                Int[]
            end
            if ep1_valid
                new_thresh = max(dot1, new_thresh)
            end
            if ep2_valid
                new_thresh = max(dot2, new_thresh)
            end
            ep1_valid &= isempty(cross1)
            ep2_valid &= isempty(cross2)
            if !(ep1_valid || ep2_valid)
                union!(remset2, cross1)
                union!(remset2, cross2)
                continue
            end
            p_new = ep1_valid ? ep1 : ep2

            push!(mesh, (p1, p2, p_new))
            # delete other edges
            eg1_exists = false
            eg2_exists = false
            try
                eidx1 = ws.edges[(p1, p_new)]
                delete!(ws, eidx1)
                eg1_exists = true
            end
            try
                eidx2 = ws.edges[(p2, p_new)]
                delete!(ws, eidx2)
                eg2_exists = true
            end
            dir_p1, dir_p2 = calc_directions(p1, p2, p_new)
            if !eg1_exists
                eidx_new1 = push!(ws, (p1, p_new), dir_p1)
                add_edge_section!(ws, eidx_new1, section)
            end
            if !eg2_exists
                eidx_new2 = push!(ws, (p2, p_new), dir_p2)
                add_edge_section!(ws, eidx_new2, section)
            end
            add_section!(ws, p_new, section)
            return true
        end
        remset, remset2 = remset2, remset
        empty!(remset2)
        threshold = new_thresh
    end
    # do_abort[] = true
    return true
    @assert false
    return false
end

function handle_next_boundary(model, ws, pset, mesh, p1, p2, p_new, step,
                              sec_p1, sec_p2, section, next_sec)
    pmid = (p1 + p2) / 2
    diff_vec = p_new - pmid
    step_len = abs(diff_vec)
    @assert step_len > 0
    if step_len < step * √(2)
        # if the boundary is closed enough, simply use that point
        push!(mesh, (p1, p2, p_new))
        dir_p1, dir_p2 = calc_directions(p1, p2, p_new)
        sec_new1 = sec_p1 == next_sec ? next_sec : section
        sec_new2 = sec_p2 == next_sec ? next_sec : section
        eidx_new1 = push!(ws, (p_new, p1), dir_p1)
        eidx_new2 = push!(ws, (p_new, p2), dir_p2)
        add_section!(ws, p_new, next_sec)
        add_edge_section!(ws, eidx_new1, sec_new1)
        add_edge_section!(ws, eidx_new2, sec_new2)
    else
        p_new, next_sec = get_next_point(model, pmid,
                                         (p_new - pmid) / 2,
                                         section, true)
        @assert next_sec == section
        push!(mesh, (p1, p2, p_new))
        dir_p1, dir_p2 = calc_directions(p1, p2, p_new)
        eidx_new1 = push!(ws, (p_new, p1), dir_p1)
        eidx_new2 = push!(ws, (p_new, p2), dir_p2)
        add_section!(ws, p_new, section)
        add_edge_section!(ws, eidx_new1, section)
        add_edge_section!(ws, eidx_new2, section)
    end
end

function handle_next_empty(model, ws, pset, mesh, p1, p2, p_new, section)
    pmid = (p1 + p2) / 2
    p_new, next_sec = get_next_point(model, pmid, (p_new - pmid) / 2,
                                     section, true)
    @assert next_sec == section
    push!(mesh, (p1, p2, p_new))
    dir_p1, dir_p2 = calc_directions(p1, p2, p_new)
    eidx_new1 = push!(ws, (p_new, p1), dir_p1)
    eidx_new2 = push!(ws, (p_new, p2), dir_p2)
    add_section!(ws, p_new, section)
    add_edge_section!(ws, eidx_new1, section)
    add_edge_section!(ws, eidx_new2, section)
end

function meshgen{V}(model::Abstract2D{V})
    pset = PointSet{V}()
    mesh = TileSet{3}(pset)
    ws = WorkSet2D{V}(pset)
    init_mesh(model, mesh, ws)

    counter = 0
    while !isempty(ws)
        (eidx, (p1, p2), dir) = first(ws)
        (sec_p1, sec_p2, section) = get_section(ws, p1, p2)
        delete!(ws, eidx)
        counter += 1
        # debug only
        # counter >= 360 && break
        # do_abort[] && break
        section == -1 && continue
        step1 = get_step_size(model, p1, section)
        step2 = get_step_size(model, p2, section)
        @assert step1 > 0
        @assert step2 > 0
        dist_p12 = abs(p1 - p2)
        step = (step1 + step2) / 2
        if step < dist_p12 / 3
            step = oftype(step, dist_p12 / 3)
        end

        pmid = (p1 + p2) / 2
        p_new, sec_new = get_next_point(model, pmid,
                                        step * oftype(step, √(2)) * dir,
                                        section, true)
        if sec_new == -1 && (p_new - pmid) * dir <= 0.1 * step
            continue
        end

        # 1. First check if the neighbering edges are good candidates, i.e.
        #
        #     1. on the correct side of the line
        #     2. small enough angle
        #     3. connecting them do not cross with other lines
        check_neighbers(model, ws, pset, mesh, p1, p2, dir, section) && continue

        # 2. Check if
        check_free_2x(model, ws, pset, mesh, p1, p2, step,
                      dir, section) && continue

        # 3. If none of the neighbering edges are good and there might be other
        #     points nearby check more carefully if there's any unterminated
        #     edges nearby. If yes, find a point among them.
        check_in_circle(model, ws, pset, mesh, p1, p2, step,
                        dir, section) && continue

        # 4. Now we know there isn't any existing points nearby, we can
        #     treat it as empty space. We first need to check if we will hit
        #     a boundary soon. Go twice as far (or some distance we know is
        #     not enough to hit any other existing edges)
        #
        #     1. If we hit an edge, truncate our step and use the point on the
        #         edge.
        if sec_new != section
            if (p_new - pmid) * dir > 0.2 * step
                handle_next_boundary(model, ws, pset, mesh, p1, p2, p_new, step,
                                     sec_p1, sec_p2, section, sec_new)
                continue
            elseif sec_new == -1
                continue
            end
            p_new, section = get_next_point(model, pmid,
                                            step * oftype(step, √(2)) * dir,
                                            section, false)
        end
        #     2. Otherwise, use the normal size
        handle_next_empty(model, ws, pset, mesh, p1, p2, p_new, section)
    end

    optimize_mesh(model, mesh, ws.edges)
    return mesh, ws.edges
end

function reconnect_safe(ep1, ep2, mp1, mp2)
    # Decide whether it's safe to swap the diagnal.
    # It's safe if the angle between the current diagnal and all the edges
    # are no larger than 90deg
    (ep1 - ep2) * (ep1 - mp1) < 0 && return false
    (ep1 - ep2) * (ep1 - mp2) < 0 && return false
    (ep2 - ep1) * (ep2 - mp1) < 0 && return false
    (ep2 - ep1) * (ep2 - mp2) < 0 && return false
    return true
end

function trig_area(p1, p2, p3)
    # Heron's Formula =)
    a = abs(p1 - p2)
    b = abs(p2 - p3)
    c = abs(p3 - p1)
    s = (a + b + c) / 2
    √(s * (s - a) * (s - b) * (s - c))
end

function reconnect_diag(midx1, midx2, ep1, ep2, mesh, edges, pset)
    mpidx1 = mesh[midx1]
    mpidx2 = mesh[midx2]

    mp1 = ep1
    mp2 = ep1

    # Find the other point in the mesh
    for r in pset[mpidx1]
        if r != ep1 && r != ep2
            mp1 = r
            break
        end
    end
    for r in pset[mpidx2]
        if r != ep1 && r != ep2
            mp2 = r
            break
        end
    end

    @assert mp1 != ep1
    @assert mp2 != ep1

    reconnect_safe(ep1, ep2, mp1, mp2) || @goto no_change

    a1 = trig_area(ep1, ep2, mp1)
    a2 = trig_area(ep1, ep2, mp2)

    a1_new = trig_area(mp1, mp2, ep1)
    a2_new = trig_area(mp1, mp2, ep2)

    a1_new * a2_new > (a1 * a2 * 1.3) && @goto do_swap
    a1 * a2 > (a1_new * a2_new * 1.3) && @goto no_change

    diag_l = abs(ep1 - ep2)
    diag_l_new = abs(mp1 - mp2)

    diag_l_new < diag_l * 0.99 && @goto do_swap

    @label no_change
    return mpidx1, mpidx2, ep1, ep2, false
    @label do_swap
    mp1_idx = pset[mp1]
    mp2_idx = pset[mp2]
    ep1_idx = pset[ep1]
    ep2_idx = pset[ep2]
    return ((ep1_idx, mp1_idx, mp2_idx), (ep2_idx, mp1_idx, mp2_idx),
            mp1, mp2, true)
end

function e2m_add(e2m, midx, mps, edges, pset, removes=())
    r1, r2, r3 = pset[mps]
    eidx1 = edges[(r1, r2)]
    eidx2 = edges[(r1, r3)]
    eidx3 = edges[(r2, r3)]
    push!(get!(e2m, eidx1, Set{Int}()), midx)
    push!(get!(e2m, eidx2, Set{Int}()), midx)
    push!(get!(e2m, eidx3, Set{Int}()), midx)
    for old_mid in removes
        delete!(e2m[eidx1], old_mid)
        delete!(e2m[eidx2], old_mid)
        delete!(e2m[eidx3], old_mid)
    end
end

function optimize_mesh(model, mesh, edges)
    pset = mesh.pts
    e2m = Dict{Int,Set{Int}}()
    for (midx, mps) in mesh.tiles
        e2m_add(e2m, midx, mps, edges, pset)
    end
    for i in 1:100
        eidxs = collect(keys(edges.tiles))
        opt_count = 0
        for eidx in eidxs
            eidx in keys(e2m) || continue
            ms = e2m[eidx]
            length(ms) == 2 || continue

            edg = edges.tiles[eidx]
            ep1, ep2 = pset[edg]
            midx1, midx2 = ms

            m1_new, m2_new, ep1, ep2, changed = reconnect_diag(midx1, midx2,
                                                               ep1, ep2,
                                                               mesh, edges, pset)

            if changed
                opt_count += 1
                delete!(mesh, midx1)
                delete!(mesh, midx2)
                delete!(edges, eidx)

                midx1_new = push!(mesh, m1_new)
                midx2_new = push!(mesh, m2_new)
                eidx = push!(edges, (ep1, ep2))
                # Fix e2m
                e2m_add(e2m, midx1_new, m1_new, edges, pset, (midx1, midx2))
                e2m_add(e2m, midx2_new, m2_new, edges, pset, (midx1, midx2))
            end
        end
        opt_count == 0 && break
        # println("optimized $opt_count edges")
    end
    mesh, edges
end

end
