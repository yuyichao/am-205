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
@generated Base.getindex{N,V}(s::PointSet{V}, vs::NTuple{N,V}) =
    Expr(:tuple, [:(s[vs[$i]]) for i in 1:N]...)
@generated Base.getindex{N,V}(s::PointSet{V}, idxs::NTuple{N,Int}) =
    Expr(:tuple, [:(s[idxs[$i]]) for i in 1:N]...)
Base.in{V<:Vec}(s::PointSet{V}, i::Int) = 1 <= i <= length(s.ary)
Base.in{V<:Vec}(s::PointSet{V}, v::V) = v in keys(s.idxs)

tuple_sort(t::Tuple{Int}) = t
tuple_sort(t::NTuple{2,Int}) = t[1] > t[2] ? (t[2], t[1]) : t
function tuple_sort(t::NTuple{3,Int})
    if t[1] <= t[2]
        if t[2] <= t[3]
            t
        elseif t[1] <= t[3]
            (t[1], t[3], t[2])
        else
            (t[3], t[1], t[2])
        end
    elseif t[1] <= t[3]
        t[2], t[1], t[3]
    elseif t[2] <= t[3]
        (t[2], t[3], t[1])
    else
        (t[3], t[2], t[1])
    end
end
function tuple_sort(t::NTuple{4,Int})
    t1 = tuple_sort((t[1], t[2]))
    t2 = tuple_sort((t[3], t[4]))
    if t1[1] <= t2[1]
        if t1[2] <= t2[1]
            (t1[1], t1[2], t2[1], t2[2])
        elseif t1[2] <= t2[2]
            (t1[1], t2[1], t1[2], t2[2])
        else
            (t1[1], t2[1], t2[2], t1[2])
        end
    elseif t2[2] <= t1[1]
        (t2[1], t2[2], t1[1], t1[2])
    elseif t2[2] <= t1[2]
        (t2[1], t1[1], t2[2], t1[2])
    else
        (t2[1], t1[1], t1[2], t2[2])
    end
end

# Slow generic version
function tuple_sort{N}(t::NTuple{N,Int})
    ary = Int[t...]
    sort!(ary)
    (ary...)::NTuple{N,Int}
end

immutable TileSet{Ndim,V<:Vec}
    pts::PointSet{V}
    cnt::Base.RefValue{Int}
    tiles::Dict{Int,NTuple{Ndim,Int}}
    idxs::Dict{NTuple{Ndim,Int},Int}
    TileSet(pts::PointSet{V}) = new(pts, Ref(0), Dict{Int,NTuple{Ndim,Int}}(),
                                    Dict{NTuple{Ndim,Int},Int}())
end
Base.call{Ndim,V}(::Type{TileSet{Ndim}}, pts::PointSet{V}) =
    TileSet{Ndim,V}(pts)
function Base.empty!(s::TileSet)
    empty!(s.tiles)
    empty!(s.idxs)
end

@generated function Base.getindex{Ndim,V}(s::TileSet{Ndim,V},
                                          tile::NTuple{Ndim,V})
    quote
        tile_idxs = ($([:(s.pts[tile[$i]]) for i in 1:Ndim]...),)
        s[tile_idxs]
    end
end
Base.getindex{Ndim,V}(s::TileSet{Ndim,V}, tile_idxs::NTuple{Ndim,Int}) =
    s.idxs[tuple_sort(tile_idxs)]
Base.getindex{Ndim,V}(s::TileSet{Ndim,V}, i::Int) = s.tiles[i]

Base.in{Ndim,V<:Vec}(s::TileSet{Ndim,V}, tile::Tuple{}) = error()
@generated function Base.in{Ndim,V<:Vec}(s::TileSet{Ndim,V},
                                         tile::NTuple{Ndim,V})
    quote
        tile_idxs = ($([:(s.pts[tile[$i]]) for i in 1:Ndim]...),)
        tile_idxs in s
    end
end
Base.in{Ndim,V<:Vec}(s::TileSet{Ndim,V}, tile_idxs::NTuple{Ndim,Int}) =
    tuple_sort(tile_idxs) in keys(s.idxs)
Base.in{Ndim,V<:Vec}(s::TileSet{Ndim,V}, i::Int) = i in keys(s.tiles)

@generated function Base.delete!{Ndim,V}(s::TileSet{Ndim,V},
                                         tile::NTuple{Ndim,V})
    quote
        tile_idxs = ($([:(s.pts[tile[$i]]) for i in 1:Ndim]...),)
        delete!(s, tile_idxs)
    end
end
function Base.delete!{Ndim,V}(s::TileSet{Ndim,V}, tile_idxs::NTuple{Ndim,Int})
    tile_idxs = tuple_sort(tile_idxs)
    tile_idxs in keys(s.idxs) || return
    idx = pop!(s.idxs, tile_idxs)
    delete!(s.tiles, idx)
    return
end
function Base.delete!{Ndim,V}(s::TileSet{Ndim,V}, i::Int)
    i in keys(s.tiles) || return
    tile_idx = pop!(s.tiles, i)
    delete!(s.idxs, tile_idx)
    return
end

@generated function Base.push!{Ndim,V}(s::TileSet{Ndim,V},
                                       tile::NTuple{Ndim,V})
    quote
        tile_idxs = tuple_sort(($([:(push!(s.pts, tile[$i]))
                                   for i in 1:Ndim]...),))
        push!(s, tile_idxs)
    end
end
@generated function Base.push!{Ndim,V}(s::TileSet{Ndim,V},
                                       tile_idxs::NTuple{Ndim,Int})
    quote
        tile_idxs = tuple_sort(tile_idxs)
        tile_idxs in keys(s.idxs) && return s.idxs[tile_idxs]
        idx = s.cnt[] + 1
        s.cnt[] = idx
        s.tiles[idx] = tile_idxs
        s.idxs[tile_idxs] = idx
        idx
    end
end


end
