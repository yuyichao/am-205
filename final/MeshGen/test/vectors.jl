#!/usr/bin/julia -f

module TestPoints

using Base.Test
import MeshGen: Vec

let
    @test isa(Vec(0.0, 1.0), Vec{2,Float64})
    @test isa(Vec((0f0, 1f0)), Vec{2,Float32})

    e_x = Vec(1.0, 0.0, 0.0)
    e_y = Vec(0.0, 1.0, 0.0)
    e_z = Vec(0.0, 0.0, 1.0)

    # == and isequal
    @test e_x == e_x
    @test e_y == e_y
    @test e_z == e_z
    @test e_x === e_x
    @test e_y === e_y
    @test e_z === e_z
    @test isequal(e_x, e_x)
    @test isequal(e_y, e_y)
    @test isequal(e_z, e_z)

    @test e_x != e_y
    @test e_y != e_z
    @test e_x != e_z

    # ==, isqual and NaN, ±0.0
    @test Vec(NaN, 0.0, 0.0) != Vec(NaN, 0.0, 0.0)
    @test !(Vec(NaN, 0.0, 0.0) == Vec(NaN, 0.0, 0.0))
    @test isequal(Vec(NaN, 0.0, 0.0), Vec(NaN, 0.0, 0.0))
    @test Vec(1.0, -0.0, 0.0) == Vec(1.0, 0.0, -0.0)
    @test Vec(1.0, -0.0, 0.0) !== Vec(1.0, 0.0, -0.0)

    # Sum and negation
    @test e_x + e_y == Vec(1.0, 1.0, 0.0) == e_y + e_x
    @test e_x + e_z == Vec(1.0, 0.0, 1.0) == e_z + e_x
    @test e_z + e_y == Vec(0.0, 1.0, 1.0) == e_y + e_z

    @test e_x + e_y + e_z == Vec(1.0, 1.0, 1.0)

    @test +e_x == e_x == -(-e_x)
    @test +e_y == e_y == -(-e_y)
    @test +e_z == e_z == -(-e_z)

    @test -e_x == Vec(-1.0, 0.0, 0.0)
    @test -e_y == Vec(0.0, -1.0, 0.0)
    @test -e_z == Vec(0.0, 0.0, -1.0)

    # Difference
    @test e_x - e_y == Vec(1.0, -1.0, 0.0)
    @test e_y - e_z == Vec(0.0, 1.0, -1.0)
    @test e_z - e_x == Vec(-1.0, 0.0, 1.0)

    # Scalar product
    @test -1.0 * e_x + 2 * e_y + e_z * 3.0 == Vec(-1.0, 2.0, 3.0)

    # Scalar product
    @test e_x / -1 + 0.5 \ e_y + e_z / (1 / 3) == Vec(-1.0, 2.0, 3.0)

    @test e_x * e_x == 1.0
    @test e_y * e_y == 1.0
    @test e_z * e_z == 1.0

    @test 2.0 * e_x * e_x == 2.0

    @test abs2(Vec{3,Complex128}(1.0im, 1.0, 0)) == 2.0
    @test abs(Vec{3,Complex128}(1.0im, 1.0, 0)) == sqrt(2.0)
end

end
