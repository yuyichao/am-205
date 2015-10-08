#!/usr/bin/julia -f

import Base: *, +, /, \, -, <, >, <=, >=

# Implementing BinNum
immutable BinNum <: Integer
    v::Bool
end

# Many of the following functions are only needed in order to use the Base
# solvers
Base.show(io::IO, n::BinNum) = Base.write(io, n.v ? "1" : "0")
Base.promote_rule{T<:Integer}(::Type{BinNum}, ::Type{T}) = BinNum

Base.convert(::Type{BinNum}, v::Bool) = BinNum(v)
Base.convert(::Type{BinNum}, v::BinNum) = v
Base.convert{T<:Integer}(::Type{BinNum}, v::T) = BinNum(v % 2 != 0)
Base.convert(::Type{Bool}, v::BinNum) = v.v
Base.convert{T<:Integer}(::Type{T}, v::BinNum) = convert(T, v.v)

Base.zero(::Type{BinNum}) = BinNum(false)
Base.one(::Type{BinNum}) = BinNum(true)
Base.zero(::BinNum) = BinNum(false)
Base.one(::BinNum) = BinNum(true)
Base.real(n::BinNum) = n
Base.abs(n::BinNum) = n
Base.isless(n1::BinNum, n2::BinNum) = isless(n1.v, n2.v)
Base.inv(n::BinNum) = BinNum(true) / n
Base.typemax(::Type{BinNum}) = BinNum(true)
Base.typemin(::Type{BinNum}) = BinNum(false)

@inline *(x::BinNum, y::BinNum) = BinNum(x.v && y.v)
@inline +(x::BinNum, y::BinNum) = BinNum(x.v $ y.v)
@inline -(x::BinNum, y::BinNum) = BinNum(x.v $ y.v)
@inline /(x::BinNum, y::BinNum) = BinNum(x.v ÷ y.v)
@inline \(x::BinNum, y::BinNum) = BinNum(y.v ÷ x.v)
@inline <(x::BinNum, y::BinNum) = x.v < y.v
@inline <=(x::BinNum, y::BinNum) = x.v <= y.v
@inline >(x::BinNum, y::BinNum) = x.v > y.v
@inline >=(x::BinNum, y::BinNum) = x.v >= y.v

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

function Base.lu(A::Matrix{BinNum})
    m, n = size(A)
    L = eye(BinNum, m)
    U = copy(A)
    Plist = [1:m;]
    l = min(m, n)
    @inbounds for j in 1:(l - 1)
        # Find a `1` in this column
        k = -1
        for i in j:m
            if U[i, j].v
                k = i
                break
            end
        end
        # Give up if there isn't one
        k == -1 && continue
        # Only do the copying if it is necessary
        if k != j
            # Update permutation list
            Plist[j], Plist[k] = Plist[k], Plist[j]
            # Update U
            for i in j:n
                U[j, i], U[k, i] = U[k, i], U[j, i]
            end
            # Update L
            for i in 1:(j - 1)
                L[j, i], L[k, i] = L[k, i], L[j, i]
            end
        end
        # Tweak the algorithm so that the inner loop is along the stride-1
        # direction
        for k in (j + 1):n
            if U[j, k].v
                # SIMD probably doesn't work either...
                @simd for i in (j + 1):m
                    U[i, k] -= U[i, j]
                end
            end
        end
        @simd for i in (j + 1):m
            L[i, j] = U[i, j]
            U[i, j] = false
        end
    end
    L, U, Plist
end

function lu_perm(A::Matrix{BinNum})
    m, n = size(A)
    L, U, Plist = lu(A)
    P = zeros(BinNum, m, m)
    @inbounds for i in 1:m
        P[i, Plist[i]] = true
    end
    L, U, P
end

# Debugging only =)
function lu_perm_base(A::Matrix{BinNum})
    m, n = size(A)
    # add the optional parameter so that we dispatch to the generic
    # Base version
    L, U, Plist = lu(A, Val{true})
    P = zeros(BinNum, m, m)
    @inbounds for i in 1:m
        P[i, Plist[i]] = true
    end
    L, U, P
end

# A = BinNum[1 0 0 1
#            1 1 0 0
#            0 1 0 0
#            0 0 1 1]

# L1, U1, P1 = lu_perm(A)
# L2, U2, P2 = lu_perm_base(A)

# println(L1 - L2)
# println(U1 - U2)
# println(P1 - P2)

function Base.tryparse_internal(::Type{BinNum}, s::AbstractString, startpos::Int,
                                endpos::Int, base::Int, raise::Bool)
    int_res = Base.tryparse_internal(Int, s, startpos, endpos, base, raise)
    !isnull(int_res) && return Nullable(BinNum(get(int_res)))
    return Nullable{BinNum}()
end

# println(readdlm("q2_large/a.txt", ' ', BinNum))
# Test for JuliaLang/julia#13483
# println(readdlm("q2_large/a.txt", ' ', BigInt))

# Which is basically `Base.(:\)` ...
function solve_lu(A, b)
    L, U, P = lu_perm(A)
    rsolve!(U, fsolve!(L, P * b))
end
