#!/usr/bin/julia -f

function generate_g!(G)
    n = size(G, 1)
    @inbounds for i in 1:(n - 1)
        # Auto vectorization already works here...
        @simd for j in 1:n
            G[j, i] = ifelse(j < i, 0, ifelse(j == i, 1, -1))
        end
    end
    # Auto vectorization already works here...
    @inbounds @simd for j in 1:n
        G[j, n] = 1
    end
    G
end

generate_g_c!(G) =
    ccall((:generate_g, "./libgen_g.so"),
          Void, (Ptr{Void}, Csize_t), G, size(G, 1))

generate_g_cilk!(G) =
    ccall((:generate_g_cilk, "./libgen_g.so"),
          Void, (Ptr{Void}, Csize_t), G, size(G, 1))

generate_g_omp!(G) =
    ccall((:generate_g_omp, "./libgen_g.so"),
          Void, (Ptr{Void}, Csize_t), G, size(G, 1))

fill_memset(G) =
    ccall(:memset, Ptr{Void}, (Ptr{Void}, Cint, Csize_t), G, 0, sizeof(G))

fill_cilk(G) =
    ccall((:fill_cilk, "./libgen_g.so"),
          Void, (Ptr{Void}, Csize_t), G, length(G))

fill_omp(G) =
    ccall((:fill_omp, "./libgen_g.so"),
          Void, (Ptr{Void}, Csize_t), G, length(G))

function generate_g(n)
    G = Array{Float64}(n, n)
    generate_g!(G)
end
