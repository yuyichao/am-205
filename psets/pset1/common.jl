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

# Periodic cubic splines

# * c_0
#    c_0   = 3x^2 - 2x^3; [0,  1]
#    c_0'  = 6x - 6x^2;   [0,  0]
#    c_0'' = 6 - 12x;     [6, -6]
# * c_1
#    c_1   = -x^2 + x^3;  [0,  0]
#    c_1'  = -2x + 3x^2;  [0,  1]
#    c_1'' = -2 + 6x;     [-2, 4]
# * c_2
#    c_2   = x^3 - 2x^2 + x; [0,  0]
#    c_2'  = 3x^2 - 4x + 1;  [1,  0]
#    c_2'' = 6x - 4;         [-4, 2]
# * c_3
#    c_3   = 2x^3 - 3x^2 + 1; [1,  0]
#    c_3'  = 6x^2 - 6x;       [0,  0]
#    c_3'' = 12x - 6;         [-6, 6]

export PeriodicSpline

immutable PeriodicSpline{T<:Number}
    dx::T
    vs::Vector{T} # Values
    ps::Vector{T} # Derivatives
    function PeriodicSpline(dx::T, vs::Vector{T})
        len = length(vs)
        if len < 2
            ArgumentError("Cannot interpolate one data point")
        end
        # First calculate the negative discontinuity caused by `vs`
        # i.e. s''(xᵢ₋) - s''(xᵢ₊) for each i if we only use c_0 and c_3
        dd2 = Vector{T}(len)
        @inbounds dd2[1] = (vs[len] - vs[2]) * 6
        @inbounds dd2[len] = (vs[len - 1] - vs[1]) * 6
        @inbounds @simd for i in 2:(len - 1)
            dd2[i] = (vs[i - 1] - vs[i + 1]) * 6
        end
        # Then calculate the contribution of the second derivative
        # discontinuity caused by each `ps`
        mat = zeros(T, len, len)
        @inbounds for i in 1:len
            mat[((i - 2 + len) % len) + 1, i] = -2
            mat[i, i] = -8
            mat[(i % len) + 1, i] = -2
        end
        ps = mat \ dd2
        new(dx, vs, ps)
    end
end

function call{T1,T2}(::Type{PeriodicSpline}, _dx::T1, vs::Vector{T2})
    dx = float(_dx)
    T = promote_type(typeof(dx), T2)
    PeriodicSpline{T}(T(dx), T[v for v in vs])
end

function call(sp::PeriodicSpline, x)
    len = length(sp.vs)
    scale_x = x / sp.dx

    num = trunc(Int, scale_x)
    x_off = scale_x - num
    idx = (num + 1) % len
    idx_lo = idx == 0 ? len : idx
    idx_hi = idx + 1

    c_0 = @evalpoly x_off 0 0 3 -2
    c_1 = @evalpoly x_off 0 0 -1 1
    c_2 = @evalpoly x_off 0 1 -2 1
    c_3 = @evalpoly x_off 1 0 -3 2

    (c_3 * sp.vs[idx_lo] + c_0 * sp.vs[idx_hi] +
     c_2 * sp.ps[idx_lo] + c_1 * sp.ps[idx_hi])
end

call(sp::PeriodicSpline, ary::AbstractVector) = [sp(x) for x in ary]

export enclosed_area

function enclosed_area(xs, ys)
    len = length(xs)
    if length(ys) != len
        throw(ArgumentError("x and y should have the same length"))
    elseif len < 2
        throw(ArgumentError("Cannot calculate enclosed area of a single point"))
    end
    T = promote_type(eltype(xs), eltype(ys))
    s::T = 0
    prev_x = xs[1]
    prev_y = ys[1]
    @inbounds for i in 2:len
        cur_x = xs[i]
        cur_y = ys[i]
        dx = cur_x - prev_x
        y = (cur_y + prev_y) / 2
        prev_x = cur_x
        prev_y = cur_y
        s += dx * y
    end
    s
end

end
