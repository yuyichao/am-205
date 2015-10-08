#!/usr/bin/julia -f

macro time_func(ex::Expr, rep)
    ex.head == :call || throw(ArgumentError("Expression has to be a function call"))
    f = ex.args[1]

    # Only works in global scope
    args = eval(current_module(), Expr(:tuple, ex.args[2:end]...))::Tuple
    argnames = Symbol[gensym("args") for i in 1:length(args)]
    types = map(typeof, args)
    quote
        function timing_wrapper()
            println($f, $types)
            $f($(args...))
            gc()
            n = 100
            times = Array{Float64}(n)
            _rep::Int = $rep
            for i in 1:n
                gc()
                t = @elapsed for j in 1:_rep
                    $f($(args...))
                end
                @inbounds times[i] = t / _rep
            end
            gc()
            mean_t = mean(times)
            mean_t2 = mean(times.^2)
            std_t = sqrt((mean_t2 - mean_t^2) / (n - 1))
            println("$mean_t ± $std_t")
            mean_t, std_t
        end
        timing_wrapper()
    end
end

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

G = zeros(1000, 1000)

@time_func generate_g!(G) 100
@time_func generate_g_c!(G) 100
@time_func generate_g_cilk!(G) 100
@time_func generate_g_omp!(G) 100
@time_func fill!(G, 0) 100
@time_func fill_memset(G) 100
@time_func fill_cilk(G) 100
@time_func fill_omp(G) 100
