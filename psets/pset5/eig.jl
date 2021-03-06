#!/usr/bin/julia -f

using PyPlot

function eigs_ode{T}(V, xs::Range{T}, n=5)
    len = length(xs)
    dx = step(xs)
    dx⁻² = 1 / dx^2
    d = Vector{T}(len)
    @inbounds for i in 1:len
        x = xs[i]
        v = V(x)
        d[i] = v + 2 * dx⁻²
    end
    du = fill(-dx⁻², len - 1)
    # Fill in a big enough value so that the edge will be zero for sure =)
    d[1] = d[len] = norm(d)
    du[1] = du[len - 1] = 0
    H = SymTridiagonal(d, du)
    eigen = eigfact(H, 1:n)
    eigen[:values], eigen[:vectors]
end

function plot_eigs{T}(V, xs::Range{T}, Ds, Vecs)
    xsize = xs[end] - xs[1]
    for i in 1:length(Ds)
        vec = slice(Vecs, (:, i))
        plot(xs, Ds[i] + vec * (xsize / 3),
             label=@sprintf("\$E_%d=%.6f\$", i, Ds[i]))
    end
    for i in 1:length(Ds)
        axhline(Ds[i], color="k")
    end
    len = length(xs)
    vs = Vector{T}(len)
    @inbounds for i in 1:len
        vs[i] = V(xs[i])
    end
    plot(xs, vs, "k", linewidth=2)
    ylim(0, Ds[end] + 2)
    xlim(xs[1], xs[end])
    grid()
    legend(loc="center left", bbox_to_anchor=(1, 0.5))
end

function int_region{T}(Vecs, xs::Range{T}, xmin, xmax)
    idxmin = max(1, round(Int, (xmin - first(xs)) / step(xs)))
    idxmax = min(length(xs), round(Int, (xmax - first(xs)) / step(xs)))
    if (idxmax - idxmin) % 2 != 1
        idxmax -= 1
    end

    ints = Vector{T}(size(Vecs, 2))

    @inbounds for i in 1:length(ints)
        int_v = Vecs[idxmin, i] + Vecs[idxmax, i]
        @simd for j in (idxmin + 1):(idxmax - 1)
            int_v = muladd(abs2(Vecs[j, i]), ifelse(j & 0x1 == 0, 4, 2), int_v)
        end
        ints[i] = int_v / 3
    end
    ints
end
