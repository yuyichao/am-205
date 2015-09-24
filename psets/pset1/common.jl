#!/usr/bin/julia -f

using PyPlot

module PSet1Common

## Interpolations

"""
An interpolation basis.

`call(base::InterpolateBase, i, x)` is used to calculate the value of
the `i`th base at point `x`.
"""
abstract InterpolateBase

"""
Interpolation result. Contains the basis and the coefficient of each base
"""
immutable InterpolationResult{B<:InterpolateBase,C<:AbstractVector}
    base::B
    coeff::C
end

"""
Calculate the value of the interpolation at `x`
"""
function call{T}(res::InterpolationResult, x::T)
    v::typeof(float(zero(T))) = 0
    len = length(res.coeff)
    @inbounds for i in 1:len
        v += res.coeff[i] * res.base(i, x)
    end
    v
end

call(res::InterpolationResult, ary::AbstractVector) = [res(x) for x in ary]

function call(res::InterpolationResult, ary::AbstractVector, out::AbstractVector)
    len = length(out)
    if length(ary) != len
        throw(ArgumentError("input and output must have the same length"))
    end
    @inbounds for i in 1:len
        out[i] = res(ary[i])
    end
end

export solve_interpolation

"""
Find the interpolation through the points given by `x` and `y` using
the basis `base_gen`.
"""
function solve_interpolation(base_gen::InterpolateBase, x::AbstractVector,
                             y::AbstractVector)
    len = length(x)
    if length(y) != len
        throw(ArgumentError("x and y must have the same length"))
    end
    Txy = promote_type(eltype(x), eltype(y))
    T = typeof(float(zero(Txy)))
    mat = Matrix{T}(len, len)
    @inbounds for i in 1:len
        for j in 1:len
            mat[j, i] = base_gen(i, x[j])
        end
    end
    InterpolationResult(base_gen, mat \ y)
end

export MonomialBase

"""
Simple monomial basis
"""
type MonomialBase <: InterpolateBase
end
call(::MonomialBase, i, x) = x^(i - 1)

export LagrangeBase

"""
Lagrange basis for a set of interpolation points
"""
immutable LagrangeBase{T<:AbstractFloat} <: InterpolateBase
    xs::Vector{T}
    denoms::Vector{T}
    function LagrangeBase(xs::Vector{T})
        len = length(xs)
        denoms = Vector{T}(len)
        @inbounds for i in 1:len
            denom::T = 1
            xi = xs[i]
            for j in 1:len
                j == i && continue
                denom *= xs[j] - xi
            end
            denoms[i] = denom
        end
        new(xs, denoms)
    end
end

call{T<:AbstractFloat}(::Type{LagrangeBase}, xs::Vector{T}) = LagrangeBase{T}(xs)
call(::Type{LagrangeBase}, xs::Vector) = LagrangeBase(float(xs))

function call{T}(base::LagrangeBase{T}, i, x)
    xs = base.xs
    len = length(xs)
    @inbounds xi = xs[i]
    @inbounds v = 1 / base.denoms[i]
    @inbounds for j in 1:len
        i == j && continue
        v *= (x - xs[j])
    end
    v
end

end
