#!/usr/bin/julia -f

import FastGaussQuadrature
using DualNumbers
using Optim

const pos251, weight251 = FastGaussQuadrature.gausslegendre(251)

@generated function parametrized_int_diff!{N,T}(f, xmin, xmax,
                                                params::NTuple{N,T}, out)
    # This hits the MAX_TUPLE_LENGTH limit ...............
    param_diff = [gensym("param_diff") for i in 1:N]
    int_res = [gensym("int_res") for i in 1:N]
    init = quote
        $([:($(param_diff[i]) = ($([:(params[$j]) for j in 1:(i - 1)]...),
                                 dual(params[$i], one(T)),
                                 $([:(params[$j]) for j in (i + 1):N]...)))
           for i in 1:N]...)
        $([:($(int_res[i])::Dual{T} = zero(Dual{T})) for i in 1:N]...)
    end
    body = quote
        @inbounds for j in 1:251
            p = pos251[j]
            x = xmax * (p + 1) / 2 - xmin * (p - 1) / 2
            $([:($(int_res[i]) += weight251[j] * f(x, $(param_diff[i])))
               for i in 1:N]...)
        end
    end
    res = quote
        scale = (xmax - xmin) / 2
        @inbounds begin
            $([:(out[$i] = $(int_res[i]).du * scale) for i in 1:N]...)
        end
    end
    quote
        $init
        $body
        $res

        out
    end
end

@generated function parametrized_int{N,T}(f, xmin, xmax, params::NTuple{N,T})
    init = quote
        int_res = zero(T)
    end
    body = quote
        @inbounds for j in 1:251
            p = pos251[j]
            x = xmax * (p + 1) / 2 - xmin * (p - 1) / 2
            int_res += weight251[j] * f(x, params)
        end
    end
    res = quote
        scale = (xmax - xmin) / 2
        int_res * scale
    end
    quote
        $init
        $body
        $res
    end
end

immutable KernelP2{N,T}
    R::T
    ω::T
    L::T
    ρ::T
    μ::T
end

@generated function Base.call{N,T}(kernel::KernelP2{N,T}, s, bs::Tuple)
    if N * 2 != length(bs.parameters)
        return :(error("...."))
    end
    quote
        R = kernel.R
        ω = kernel.ω
        L = kernel.L
        ρ = kernel.ρ
        μ = kernel.μ
        y = +($([:((bs[N + $i]::$(bs.parameters[N + i])) * sin(π * $i * s / R))
                 for i in 1:N]...))
        dy = π * +($([:((bs[N + $i]::$(bs.parameters[N + i])) * $i *
                         cos(π * $i * s / R)) for i in 1:N]...))
        dx = L / R + π * +($([:((bs[$i]::$(bs.parameters[i])) * $i *
                                 cos(π * $i * s / R)) for i in 1:N]...))
        μ * (√(dx^2 + dy^2) - 1)^2 - ρ * y^2 * ω^2
    end
end

immutable OptimIntKernel{N,F,T}
    f::F
    xmin::T
    xmax::T
end
call{N,F,T}(::Type{OptimIntKernel{N}}, f::F, xmin::T, xmax::T) =
    OptimIntKernel{N,F,T}(f, xmin, xmax)
@generated function call{N,F,T}(k::OptimIntKernel{N,F,T}, x)
    :(parametrized_int(k.f, k.xmin, k.xmax, ($([:(x[$i]) for i in 1:N]...),)))
end
@generated function call{N,F,T}(k::OptimIntKernel{N,F,T}, x, out)
    :(parametrized_int_diff!(k.f, k.xmin, k.xmax,
                             ($([:(x[$i]) for i in 1:N]...),), out))
end

function calc_x{N,T}(kernel::KernelP2{N,T}, s, bs::Vector{T})
    if N * 2 != length(bs)
        throw(ArgumentError("...."))
    end
    R = kernel.R
    ω = kernel.ω
    L = kernel.L
    ρ = kernel.ρ
    μ = kernel.μ
    x = L * s / R
    for i in 1:N
        x += bs[i] * sin(π * i * s / R)
    end
    x
end

function calc_y{N,T}(kernel::KernelP2{N,T}, s, bs::Vector{T})
    if N * 2 != length(bs)
        throw(ArgumentError("...."))
    end
    R = kernel.R
    ω = kernel.ω
    L = kernel.L
    ρ = kernel.ρ
    μ = kernel.μ
    y = zero(T)
    for i in 1:N
        y += bs[i + N] * sin(π * i * s / R)
    end
    y
end
