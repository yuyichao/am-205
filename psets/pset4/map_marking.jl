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
    data
end

mark_floor_map!(data, (59, 17))

imshow(data, interpolation="none")
savefig("pierce_mark_reachable.png")
show()
