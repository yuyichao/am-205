#!/usr/bin/julia -f

module P5Opt

using ForwardDiff

type LineSearchWrapper{T,F}
    f::F
    grad::Vector{T}
    x0::Vector{T}
end
LineSearchWrapper{T,F}(f::F, grad::Vector{T}, x0::Vector{T}) =
    LineSearchWrapper{T,F}(f, grad, x0)

function Base.call(w::LineSearchWrapper, _x)
    grad = w.grad
    x0 = w.x0
    x = _x[1]
    w.f(grad * x + x0)
end

function do_line_search{T}(ls::LineSearchWrapper{T})
    x = T[0]
    grad = T[0]
    while true
        # Would be nice if ForwardDiff has a method that doesn't allocate
        hess, allres = ForwardDiff.hessian(ls, x, ForwardDiff.AllResults)
        ForwardDiff.gradient!(grad, allres)
        dx = -(grad[1] / hess[1])
        x[1] += dx
        abs(dx) < 1e-8 && return x[1]
    end
end

export steepest_descent
function steepest_descent(f, r0, accum)
    line_search = LineSearchWrapper(f, similar(r0), copy(r0))
    grad = line_search.grad
    cur_r = line_search.x0
    accum(cur_r)
    while true
        gradient!(grad, f, cur_r)
        ldx = do_line_search(line_search)
        @inbounds @simd for i in eachindex(grad)
            cur_r[i] = muladd(grad[i], ldx, cur_r[i])
        end
        accum(cur_r)
        norm(grad) * abs(ldx) < 1e-8 && return cur_r
    end
end

export newton_min
function newton_min(f, r0, accum)
    cur_r = copy(r0)
    grad = similar(r0)
    accum(cur_r)
    while true
        hess, allres = ForwardDiff.hessian(f, cur_r, ForwardDiff.AllResults)
        ForwardDiff.gradient!(grad, allres)
        dr = hess \ grad
        @inbounds @simd for i in eachindex(cur_r)
            cur_r[i] -= dr[i]
        end
        accum(cur_r)
        norm(dr) < 1e-8 && return cur_r
    end
end

export bfgs_min
function bfgs_min(f, r0, accum)
    len = length(r0)

    next_grad = similar(r0)
    grad = ForwardDiff.gradient(f, r0)
    cur_r = copy(r0)
    H = eye(length(r0))
    T = similar(H)
    tmpm1 = similar(H)
    dr = similar(cur_r)

    accum(cur_r)
    @inbounds while true
        A_mul_B!(dr, H, grad)
        for i in 1:len
            cur_r[i] -= dr[i]
        end
        accum(cur_r)
        norm(dr) < 1e-8 && return cur_r
        ForwardDiff.gradient!(next_grad, f, cur_r)
        # reuse grad as y
        for i in 1:len
            grad[i] = next_grad[i] - grad[i]
        end
        ρ = -1 / dot(grad, dr)
        # T = I + ρ dr * grad'
        for j in 1:len
            yj = grad[j] * ρ
            for i in 1:len
                T[i, j] = dr[i] * yj
            end
            T[j, j] += 1
        end
        # T * H * T'
        A_mul_B!(tmpm1, T, H)
        A_mul_Bt!(H, tmpm1, T)
        # H += ρ dr * dr'
        for j in 1:len
            dr_j = dr[j] * ρ
            for i in 1:len
                H[i, j] += dr[i] * dr_j
            end
        end
        next_grad, grad = grad, next_grad
    end
end

export RecordXY
immutable RecordXY{T}
    xs::Vector{T}
    ys::Vector{T}
    RecordXY() = new(T[], T[])
end
function Base.call(record::RecordXY, r)
    push!(record.xs, r[1])
    push!(record.ys, r[2])
    nothing
end

export Rosebrock
immutable Rosebrock
end
@inline function Base.call(::Rosebrock, x, y)
    100 * (y - x^2)^2 + (1 - x)^2
end
@inline function Base.call(f::Rosebrock, r)
    f(r[1], r[2])
end

end
