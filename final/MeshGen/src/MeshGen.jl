#!/usr/bin/julia -f

module MeshGen

include("vectors.jl")
import .Vectors: Vec

include("tilesets.jl")
import .TileSets: PointSet, TileSet

include("meshes.jl")
import .Meshes: mesh

end
