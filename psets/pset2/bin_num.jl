#!/usr/bin/julia -f

import Base: *, +, /, \

immutable BinNum
    v::Bool
end

Base.show(io::IO, n::BinNum) = Base.write(io, n.v ? "1b" : "0b")

Base.convert(::Type{BinNum}, v::Bool) = BinNum(v)
Base.convert{T<:Integer}(::Type{BinNum}, v::T) = BinNum(v != 0)
Base.zero(::Type{BinNum}) = BinNum(false)
Base.one(::Type{BinNum}) = BinNum(true)
Base.zero(::BinNum) = BinNum(false)
Base.one(::BinNum) = BinNum(true)

@inline *(x::BinNum, y::BinNum) = BinNum(x.v && y.v)
@inline +(x::BinNum, y::BinNum) = BinNum(x.v $ y.v)
@inline /(x::BinNum, y::BinNum) = BinNum(x.v ÷ y.v)
@inline \(x::BinNum, y::BinNum) = BinNum(y.v ÷ x.v)

# L = BinNum[1 0 0 0
#            0 1 0 0
#            1 1 1 0
#            1 0 1 1]

# U = BinNum[1 0 1 0
#            0 1 1 1
#            0 0 1 0
#            0 0 0 1]

# println(L * U)

# Modify L, b and out in place
function fsolve!(L::Matrix{BinNum}, b::Vector{BinNum}, out=similar(b))
    len = length(b)
    if size(L) != (len, len) || length(out) != len
        throw(ArgumentError("Whatever...."))
    end
    @inbounds for i in 1:len
        if !L[i, i].v
            throw(DivideError())
        end
        if b[i].v
            # SIMD for this doesn't currently work
            # JuliaLang/julia#13104
            @simd for j in (i + 1):len
                b[j] += L[j, i]
            end
        end
        @simd for j in (i + 1):len
            L[j, i] = false
        end
    end
    b
end

# Modify U, b and out in place
function rsolve!(U::Matrix{BinNum}, b::Vector{BinNum}, out=similar(b))
    len = length(b)
    if size(U) != (len, len) || length(out) != len
        throw(ArgumentError("Whatever...."))
    end
    @inbounds for i in len:-1:1
        if !U[i, i].v
            throw(DivideError())
        end
        if b[i].v
            # SIMD for this doesn't currently work
            # JuliaLang/julia#13104
            @simd for j in 1:(i - 1)
                b[j] += U[j, i]
            end
        end
        @simd for j in 1:(i - 1)
            U[j, i] = false
        end
    end
    b
end
