#!/usr/bin/julia -f

module MeshGen

function get_section end
# step size should be no greater than the thickness of the next section
function get_step_size end
function get_edge end
function check_crossing end

abstract Abstract1D{T}

# `0` is assumed to be a point in a closed loop or the negative end of a line
function mesh_1d{T}(model::Abstract1D{T})
    mesh = T[0]
    frontier = T[0, 0]
    cur_section = get_section(model, frontier[2])
    @inbounds while true
        p = frontier[2]
        # Get the position of the next point
        step_size = get_step_size(model, p, cur_section)
        next_p = p + step_size
        # Check if we crosses some boundary
        next_section = get_section(model, next_p)
        if next_section != cur_section
            edge = get_edge(model, cur_section, next_section, p, next_p)
            # edge is a valid point, but we still need to check if it crosses
            # with the region we already meshed.
            check_crossing(model, edge, frontier) && return mesh
            push!(mesh, edge)
            # If we reaches the edge
            next_section == -1 && return mesh
            # This is a normal boundary
            frontier[2] = edge
            cur_section = next_section
            continue
        end
        check_crossing(model, next_p, frontier) && return mesh
        push!(mesh, next_p)
        frontier[2] = next_p
    end
end

end
