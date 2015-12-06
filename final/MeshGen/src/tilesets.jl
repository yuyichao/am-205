#!/usr/bin/julia -f

module TileSets

import ..Vectors: Vec

immutable PointSet{V<:Vec}
    ary::Vector{V}
    idxs::Dict{V,Int}
    PointSet() = new(V[], Dict{V,Int}())
end
Base.getindex{V}(s::PointSet{V}, i::Int) = s.ary[i]
Base.getindex{V}(s::PointSet{V}, v::V) = s.idxs[v]
function Base.push!{V}(s::PointSet{V}, v::V)
    v in keys(s.idxs) && return s.idxs[v]
    push!(s.ary, v)
    idx = length(s.ary)
    s.idxs[v] = idx
    idx
end

immutable TileSet{Ndim,V<:Vec}
    pts::PointSet{V}
    tiles::Vector{NTuple{Ndim,Int}}
    idxs::Dict{V,Int}
    TileSet(pts::PointSet{V}) = new(pts, NTuple{Ndim,Int}[])
end
Base.call{Ndim,V}(::Type{TileSet{Ndim}}, pts::PointSet{V}) =
    TileSet{Ndim,V}(pts)
@generated function Base.getindex{Ndim,V}(s::TileSet{Ndim,V},
                                          tile::NTuple{Ndim,V})
    quote
        tile_idxs = ($([:(s.pts[tile[$i]]) for i in 1:Ndim]...),)
        s[tile_idxs]
    end
end
Base.getindex{Ndim,V}(s::TileSet{Ndim,V}, tile_idxs::NTuple{Ndim,Int}) =
    s.idxs[tile_idxs]
Base.getindex{Ndim,V}(s::TileSet{Ndim,V}, i::Int) = s.tiles[i]
@generated function Base.push!{Ndim,V}(s::TileSet{Ndim,V},
                                       tile::NTuple{Ndim,V})
    quote
        tile_idxs = ($([:(push!(s.pts, tile[$i])) for i in 1:Ndim]...),)
        tile_idxs in keys(s.idxs) && return s.idxs[tile_idxs]
        push!(s.tiles, tile_idxs)
        idx = length(s.tiles)
        s.idxs[tile_idxs] = idx
        idx
    end
end


end
