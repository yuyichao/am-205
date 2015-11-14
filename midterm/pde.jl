#!/usr/bin/julia -f

# What the problem asked for is too simple so it doesn't really worth
# implementing a more generic version....
function propagate_pde_p2(u0, ν, α, nsteps)
    u_prev = copy(u0)
    u_cur = similar(u0)
    a₋₁ = -(α - 1) * ν
    a₀ = 1 - (1 - 2α) * ν
    a₁ = -α * ν
    nele = length(u0)
    @inbounds for i in 1:nsteps
        u_cur[1] = a₋₁ * u_prev[nele] + a₀ * u_prev[1] + a₁ * u_prev[2]
        u_cur[nele] = (a₋₁ * u_prev[nele - 1] +
                       a₀ * u_prev[nele] + a₁ * u_prev[1])
        @simd for i in 2:(nele - 1)
            u_cur[i] = a₋₁ * u_prev[i - 1] + a₀ * u_prev[i] + a₁ * u_prev[i + 1]
        end
        u_cur, u_prev = u_prev, u_cur
    end
    u_prev
end

get_u0(h, neles) = [1 / (2 + cospi(2 * h * i)) for i in 0:(neles - 1)]
