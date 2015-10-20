#!/usr/bin/julia -f

module Integrate

function comp_trap(func, xs)
    s = float(zero(eltype(xs)))
    len = length(xs)
    prev_x = xs[1]
    prev_y = func(prev_x)
    @inbounds for i in 2:len
        x = xs[i]
        y = func(x)
        δx = x - prev_x
        s += (y + prev_y) * δx / 2
        prev_x = x
        prev_y = y
    end
    s
end

function integrate_g3{T1,T2}(func, a::T1, b::T2)
    T = typeof(float(zero(promote_type(T1, T2))))
    sq35 = sqrt(T(3 / 5))
    x0 = (a + b) / 2
    δx = (b - a) / 2
    x₋ = x0 - δx * sq35
    x₊ = x0 + δx * sq35
    ((func(x₋) + func(x₊)) * (5 / 9) + func(x0) * (8 / 9)) * δx
end

immutable AdaptiveResult{T}
    res::T
    err::T
    count::Int
end

Base.(:+)(res1::AdaptiveResult, res2::AdaptiveResult) =
    AdaptiveResult(res1.res + res2.res, res1.err + res2.err,
                   res1.count + res2.count)

Base.show(io::IO, res::AdaptiveResult) =
    @printf(io, "Value: %.8g; Error: %.2g; Interval count: %d",
            res.res, res.err, res.count)

function _adaptive_g3(func, I0, a, b, t)
    I1 = integrate_g3(func, a, (a + b) / 2)
    I2 = integrate_g3(func, (a + b) / 2, b)
    I′ = I1 + I2
    err = abs(I0 - I′)
    if err < t * (b - a)
        return AdaptiveResult(I′, err, 2)
    end
    return (_adaptive_g3(func, I1, a, (a + b) / 2, t) +
            _adaptive_g3(func, I2, (a + b) / 2, b, t))
end

adaptive_g3(func, a, b, t=1e-6) =
    _adaptive_g3(func, integrate_g3(func, a, b), a, b, t)

function adaptive_g3_rand(func, _a, _b, t=1e-6)
    a, b = promote(float(_a), float(_b))
    T = typeof(a)
    uppers = [[a + (b - a) * x for x in sort!(rand(T, 10))]; b]
    res = adaptive_g3(func, a, uppers[1], t)
    prev_u = uppers[1]
    for i in 2:length(uppers)
        u = uppers[i]
        res += adaptive_g3(func, prev_u, u, t)
        prev_u = u
    end
    res
end

end
