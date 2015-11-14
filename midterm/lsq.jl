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
