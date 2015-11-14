#!/usr/bin/julia -f

function fit_F(x::AbstractMatrix, r::AbstractMatrix)
    @assert size(x, 2) == 2
    @assert size(r, 2) == 2
    len = size(x, 1)
    @assert size(r, 1) == len
    A = Matrix{eltype(x)}(2 * len, 6)
    b = Vector{eltype(r)}(2 * len)
    @inbounds for i in 1:len
        A[i * 2 - 1, 1] = x[i, 1]
        A[i * 2, 1] = 0
        A[i * 2 - 1, 2] = x[i, 2]
        A[i * 2, 2] = 0
        A[i * 2 - 1, 3] = 1
        A[i * 2, 3] = 0

        A[i * 2 - 1, 4] = 0
        A[i * 2, 4] = x[i, 1]
        A[i * 2 - 1, 5] = 0
        A[i * 2, 5] = x[i, 2]
        A[i * 2 - 1, 6] = 0
        A[i * 2, 6] = 1

        b[i * 2 - 1] = r[i, 1]
        b[i * 2] = r[i, 2]
    end
    fc = A \ b
    # fc = At_mul_B(A, A) \ At_mul_B(A, b)
    [fc[1] fc[2]
     fc[4] fc[5]], [fc[3], fc[6]]
end

function fit_Fc(data::AbstractMatrix)
    # Calculate F_c, plus many random stuff the problem is asked for...
    @assert size(data, 2) == 4
    len = size(data, 1)

    # get mean positions
    sum_x = sum_y = sum_r = sum_s = zero(eltype(data))
    @inbounds @simd for i in 1:len
        sum_x += data[i, 1]
        sum_y += data[i, 2]
        sum_r += data[i, 3]
        sum_s += data[i, 4]
    end
    avg_x = sum_x / len
    avg_y = sum_y / len
    avg_r = sum_r / len
    avg_s = sum_s / len

    A = similar(data)
    @inbounds @simd for i in 1:len
        A[i, 1] = data[i, 1] - avg_x
        A[i, 2] = data[i, 2] - avg_y
        A[i, 3] = data[i, 3] - avg_r
        A[i, 4] = data[i, 4] - avg_s
    end

    svd_res = svdfact(A)
    # This is not the official API and might break at any time
    U = svd_res.U
    S = svd_res.S
    Vt = svd_res.Vt

    # Not super efficient and assume S is ordered
    A′ = U * diagm([S[1], S[2], 0, 0]) * Vt

    E1 = zero(eltype(data))
    @inbounds @simd for i in 1:len
        E1 += (abs2(A′[i, 1] - A[i, 1]) + abs2(A′[i, 2] - A[i, 2]) +
               abs2(A′[i, 3] - A[i, 3]) + abs2(A′[i, 4] - A[i, 4]))
    end
    E2 = abs2(S[3]) + abs2(S[4])
    @assert E1 ≈ E2

    x′ = sub(A′, (:, 1:2))
    r′ = sub(A′, (:, 3:4))
    # Reuse the function above since why not....
    F_c, c_c = fit_F(x′, r′)
    F_c, E1
end
