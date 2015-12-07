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
end

function mesh{V}(model::Abstract1D{V})
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
                @assert step_len > 0
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

end
