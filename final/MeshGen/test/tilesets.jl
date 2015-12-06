#!/usr/bin/julia -f

module TestTileSets

using Base.Test
import MeshGen: Vec, PointSet, TileSet

let
    v1 = Vec(0.0, 1.0)
    v2 = Vec(1.0, 0.0)

    pset = PointSet{Vec{2,Float64}}()
    @test_throws KeyError pset[v1]
    @test_throws KeyError pset[v2]
    @test_throws BoundsError pset[1]

    @test push!(pset, v1) == 1
    @test push!(pset, v2) == 2
    @test push!(pset, v1) == 1
    @test push!(pset, v2) == 2

    @test pset[v1] == 1
    @test pset[v2] == 2
    @test pset[1] == v1
    @test pset[2] == v2

    @test pset[(1, 2)] == (v1, v2)
    @test pset[(v1, v2)] == (1, 2)

    v3 = Vec(0.0, 0.0)

    tset = TileSet{3}(pset)
    @test_throws KeyError tset[(v1, v2, v3)]
    @test push!(tset, (v2, v1, v3)) == 1
    @test push!(tset, (v1, v2, v3)) == 1
    @test tset[(v1, v2, v3)] == 1
    @test tset[(1, 2, 3)] == 1
    @test tset[1] == (1, 2, 3)
    @test pset[v3] == 3
    @test pset[3] == v3

    delete!(tset, (v1, v2, v3))
    @test_throws KeyError tset[(v1, v2, v3)]
end

end
