#!/usr/bin/julia -f

"""
Solve the ODE defined by `y′ = f(t, y)` using the Runge-Kutta method.
"""
function ode_rk4(f, t0, y0, h, cb)
    h2 = h / 2
    h3 = h / 3
    h6 = h / 6
    y_prev = y0
    @inbounds for i in countfrom()
        t = t0 + h * i
        k1 = f(t, y_prev)
        y = muladd(h2, k1, y_prev)
        k2 = f(t + h2, y)
        y = muladd(h2, k2, y_prev)
        k3 = f(t + h2, y)
        y = muladd(h, k3, y_prev)
        k4 = f(t + h, y)
        y = muladd(h6, k1, y_prev) + muladd(h3, k2, h3 * k3) + h6 * k4
        cb(y, i, t) || break
    end
end

type AbortNegative{T}
    t_prev::T
    y_prev::T
    t_last::T
    y_last::T
    AbortNegative() = new(NaN, NaN, NaN, NaN)
end

function call(cb::AbortNegative, y, i, t)
    if y >= 0
        cb.t_prev = t
        cb.y_prev = y
        return true
    end
    cb.t_last = t
    cb.y_last = y
    return false
end

immutable CubicSpline{T}
    a0::T
    a1::T
    a2::T
    a3::T
    function CubicSpline(y1, y2, y1′, y2′)
        #  y1:  2x^3 - 3x^2     + 1
        #  y2: -2x^3 + 3x^2
        # y1′:   x^3 - 2x^2 + x
        # y2′:   x^3 -  x^2
        a0 = y1
        a1 = y1′ / δt
        a2 = -3 * y1 + 3 * y2 - 2 * y1′ - y2′ / δt^2
        a3 = 2 * y1 - 2 * y2 + y1′ + y2′ / δt^3
        new(a0, a1, a2, a3)
    end
end

CubicSpline{T}(y1::T, y2::T, y1′::T, y2′::T) =
    CubicSpline{T}(y1, y2, y1′, y2′)

@inline call(sp::CubicSpline, x) = @evalpoly(x, sp.a0, sp.a1, sp.a2, sp.a3)

function find_ode_zero{T}(f, t0, y0::T, h)
    @assert y0 >= 0
    cb = AbortNegative{T}()
    ode_rk4(f, t0, y0, h, cb)
    @assert !isnan(cb.y_prev)
    t_prev = cb.t_prev
    y_prev = cb.y_prev
    t_last = cb.t_last
    y_last = cb.y_last
    y′_prev = f(t_prev, y_prev)
    y′_last = f(t_last, y_last)
    @assert y_prev >= 0 && y_last < 0
    sp = CubicSpline(y_prev, y_last, y_prev′, y_last′)
    # Do bisect on relative time to avoid doing floating point math on ulp level
    ɛ = oftype(h, 1e-15) / h
    x_lower = zero(t_prev)
    x_upper = one(t_last)
    x_mid = x_upper
    while x_upper - x_lower >= ɛ
        x_mid = (x_upper + x_lower) / 2
        if sp(x_mid) >= 0
            x_lower = x_mid
        else
            x_upper = x_mid
        end
    end
    x_mid * h + t_last
end
