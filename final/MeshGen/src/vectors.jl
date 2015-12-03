#!/usr/bin/julia -f

module Vectors

import Base: call, ==, isequal, +, -, *, /, \, abs, abs2

immutable Vec{N,T}
    r::NTuple{N,T}
    @generated function Vec(rs::Tuple)
        if !isleaftype(T)
            return :(error("Only leaf type is allowed. Got $T instead"))
        end
        Expr(:new, Vec{N,T}, Expr(:tuple, [:(T(rs[$i])) for i in 1:N]...))
    end
    @inline Vec(rs...) = Vec{N,T}(rs)
end
call{N,T}(::Type{Vec}, rs::NTuple{N,T}) = Vec{N,T}(rs)
call(::Type{Vec}, rs::Number...) = Vec(rs)

@generated function (==){N}(vec1::Vec{N}, vec2::Vec{N})
    expr = true
    for i in 1:N
        expr = :($expr && vec1.r[$i] == vec2.r[$i])
    end
    expr
end
@generated function isequal{N}(vec1::Vec{N}, vec2::Vec{N})
    expr = true
    for i in 1:N
        expr = :($expr && isequal(vec1.r[$i], vec2.r[$i]))
    end
    expr
end

@generated function +{N}(vec1::Vec{N}, vecs::Vec{N}...)
    len = length(vecs)
    quote
        $(Expr(:meta, :inline))
        Vec($([:(+(vec1.r[$i], $([:(vecs[$j].r[$i]) for j in 1:len]...)))
               for i in 1:N]...))
    end
end

+(vec::Vec) = vec

@generated function -{N}(vec::Vec{N})
    :(Vec($([:(-vec.r[$i]) for i in 1:N]...)))
end

@generated function -{N}(vec1::Vec{N}, vec2::Vec{N})
    :(Vec($([:(vec1.r[$i] - vec2.r[$i]) for i in 1:N]...)))
end

@generated function *{N,T<:Number}(s::T, vec::Vec{N})
    :(Vec($([:(s * vec.r[$i]) for i in 1:N]...)))
end
@generated function *{N,T<:Number}(vec::Vec{N}, s::T)
    :(Vec($([:(vec.r[$i] * s) for i in 1:N]...)))
end

@generated function /{N,T<:Number}(vec::Vec{N}, s::T)
    :(Vec($([:(vec.r[$i] / s) for i in 1:N]...)))
end
@generated function \{N,T<:Number}(s::T, vec::Vec{N})
    :(Vec($([:(s \ vec.r[$i]) for i in 1:N]...)))
end

@generated function *{N}(vec1::Vec{N}, vec2::Vec{N})
    :(+($([:(vec1.r[$i] * vec2.r[$i]) for i in 1:N]...)))
end

@generated function abs2{N}(vec::Vec{N})
    :(+($([:(abs2(vec.r[$i])) for i in 1:N]...)))
end
@inline abs(vec::Vec) = sqrt(abs2(vec))

end
