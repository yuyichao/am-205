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

end
