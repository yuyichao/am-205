#!/usr/bin/julia -f

using PyPlot

data = readdlm("pierce.txt", ' ', UInt32)

# The data is presented in a way that is not very friendly for calculation for
# a few reasons.
#
# 1. Even though `0` are always wall and therefore is sth we don't care about
#    `1` doesn't necessarily mean we care about them either since they might
#    not be connected (separated by the wall)
# 2. The `0`, `1` color also only contains information of the current ceil.
#    However, for the simulation, we'd like to know whether this is a ceil next
#    to the wall since we'd like to load values differently in such cases.
# 3. Finally (minor), the map doesn't tell us the index range we are actually
#    interested in. A large fraction of the map are walls or otherwise
#    inaccessible areas. It would be quite costy if we have to loop over all
#    those areas.
#
# Therefore, we'd like to do some pre-processing of the map such that it is
# in the format that is easiest for us to use.
#
# First, let's figure out the region we are actually intersted in.

function mark_floor_map!(data, start)
    # We start from the frontier (i.e. the points that have just been marked)
    # of the previos iteration and mark all its neighbers that are not wall.
    # We skip bounds check and assume the area is surrounded by a wall.
    frontier = Tuple{Int,Int}[]
    prev_frontier = Tuple{Int,Int}[start]
    # We use `2` to mark the points that are visited so we don't visit them
    # again
    data[start...] = 2
    @inbounds while !isempty(prev_frontier)
        # Visit all the points in the frontier
        for (x, y) in prev_frontier
            # Mark all their neighbers
            for (x′, y′) in ((x, y - 1), (x - 1, y), (x + 1, y), (x, y + 1))
                v = data[x′, y′]
                if v == 0
                    # Unvisited empty space
                    data[x′, y′] = 2
                    push!(frontier, (x′, y′))
                end
            end
        end
        # Now we have marked all the points from the previous frontier and
        # have generated a new frontier. Clear the previous one, swap them
        # and start over.
        empty!(prev_frontier)
        prev_frontier, frontier = frontier, prev_frontier
    end
end

mark_floor_map!(data, (59, 17))

# imshow(data, interpolation="none")
# savefig("pierce_mark_reachable.png")

# Next we can mark on each pixels what kind of space it is.
# We use bits flags to store whether the point, or its neighbers are valid
# grid points. In the order of self, y - 1, x - 1, x + 1, y + 1 (from low to
# high bit starting with the lowest bit)
function mark_neighbers!(data)
    # The value on each pixels are really inconsistent right now. Since we
    # aren't really going through the pixels in a simple order to mark the
    # type of their neighbers, it's easier to normalize the values first.
    @inbounds @simd for i in eachindex(data)
        v = data[i]
        data[i] = ifelse(v == 2, 1, 0)
    end
    # Now we can actually mark the neighber types
    # Again, we assume the edges are not interesting and just skip them
    # (Too lazy to write code specific for the edges that will not be
    # executed anyway....)
    @inbounds for j in 2:(size(data, 2) - 1)
        @simd for i in 2:(size(data, 1) - 1)
            # We mask off the wall while marking the neighbers so `v == 0` is
            # enough to tell if it's in the wall (instead of `v & 1 == 0`)
            v_my = ifelse(data[i, j - 1] == 0, 0x0, 0x2)
            v_mx = ifelse(data[i - 1, j] == 0, 0x0, 0x4)
            # This load is free, we write to it anyway later. Might as well use
            # it to mask off the wall completely
            v = data[i, j] == 0
            v_px = ifelse(data[i + 1, j] == 0, 0x0, 0x8)
            v_py = ifelse(data[i, j + 1] == 0, 0x0, 0x10)
            v_neighbers = (v_my | v_mx) | (v_px | v_py)
            data[i, j] = ifelse(v, 0, v_neighbers | 0x1)
        end
    end
end

mark_neighbers!(data)
# imshow(data, interpolation="none")
# savefig("pierce_mark_neighbers.png")
