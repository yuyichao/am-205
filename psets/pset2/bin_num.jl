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

*(x::BinNum, y::BinNum) = BinNum(x.v && y.v)
+(x::BinNum, y::BinNum) = BinNum(x.v $ y.v)
/(x::BinNum, y::BinNum) = BinNum(x.v ÷ y.v)
\(x::BinNum, y::BinNum) = BinNum(y.v ÷ x.v)

# L = BinNum[1 0 0 0
#            0 1 0 0
#            1 1 1 0
#            1 0 1 1]

# U = BinNum[1 0 1 0
#            0 1 1 1
#            0 0 1 0
#            0 0 0 1]

# println(L * U)
