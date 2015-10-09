#!/usr/bin/julia -f

# Q₀ᵀA = R₀
# Q′Q₀ᵀA = Q′R₀
# R₁ = Q′R₀
# Q₁ = Q₀Q′ᵀ
function givens!(A)
    m, n = size(A)
    @assert m >= n
    Q = eye(eltype(A), m)
    cur_col = Vector{eltype(A)}(m)
    cur_col2 = Vector{eltype(A)}(m)
    @inbounds for k in 1:n
        # This copy is unecessary
        @simd for i in 1:m
            cur_col[i] = A[i, k]
        end
        # This also does extra work since many of the element in Q are zero
        At_mul_B!(cur_col2, Q, cur_col)
        a1 = cur_col2[m]
        for j in m:-1:(k + 1)
            a2 = a1
            a1 = cur_col2[j - 1]
            # [ c  s] [a1]
            # [-s  c] [a2]
            c, s = if a1 > a2
                t = a2 / a1
                _c = 1 / √(1 + t^2)
                _c, _c * t
            else
                τ = a1 / a2
                _s = 1 / √(1 + τ^2)
                _s * τ, _s
            end
            # [ Qⱼ₋₁ Qⱼ  ] [ c -s]
            # [ ...  ... ] [ s  c]
            @simd for i in 1:m
                Q_1 = Q[i, j - 1]
                Q_2 = Q[i, j]
                Q_1′ = c * Q_1 + s * Q_2
                Q_2′ = -s * Q_1 + c * Q_2
                Q[i, j - 1] = Q_1′
                Q[i, j] = Q_2′
            end
            a1 = c * a1 + s * a2
        end
        cur_col2[k] = a1
        @simd for i in 1:k
            A[i, k] = cur_col2[i]
        end
        @simd for i in (k + 1):m
            A[i, k] = 0
        end
    end
    Q, A
end
