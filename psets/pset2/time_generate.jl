#!/usr/bin/julia -f

using JSON

include("generate_g.jl")

compile_typeof(a) = isa(a,Type) ? Type{a} : typeof(a)

# Use in local scope only
macro time_func(ex::Expr)
    if ex.head != :call
        throw(ArgumentError("Expression has to be a function call"))
    end
    f = ex.args[1]
    args = ex.args[2:end]
    nargs = length(args)
    arg_vars = [gensym("args") for i in 1:nargs]
    ex_str = string(ex)

    quote
        let
            # Run the benchmark for at least 10s
            min_total = 10.0
            # Run each iteration for at least 0.1s
            min_single = 0.1
            println(string("Benchmarking: ", $ex_str))
            $([:($(arg_vars[i]) = $(args[i])) for i in 1:nargs]...)
            # Precompile it, just in case
            precompile($f, ($([:(compile_typeof($(arg_vars[i])))
                               for i in 1:nargs]...),))
            gc()
            t_baseline = @elapsed $f($(arg_vars...))
            if t_baseline <= min_single
                gc()
                t_baseline = @elapsed for i in 1:10
                    $f($(arg_vars...))
                end
                t_baseline /= 10
            end
            numloop = max(round(Int, min_single / t_baseline), 1)
            t_loop = t_baseline * numloop
            # Run at least 3 times to get some statistics
            numrep = max(round(Int, min_total / t_loop), 3)
            # (Kind of) predictable starting point
            gc()
            gc()
            times = Array{Float64}(numrep)
            for i in 1:numrep
                gc()
                t = @elapsed for j in 1:numloop
                    $f($(arg_vars...))
                end
                @inbounds times[i] = t / numloop
            end
            gc()
            mean_t = mean(times)
            mean_t2 = mean(times.^2)
            std_t = sqrt((mean_t2 - mean_t^2) / (numrep - 1))
            println("    $numloop loops, repeat $numrep times")
            println("    t = $mean_t ± $std_t")
            mean_t, std_t
        end
    end
end

function time_generate(n)
    G = zeros(n, n)

    t_alloc = @time_func generate_g(n)
    t_gen = @time_func generate_g!(G)
    t_gen_c = @time_func generate_g_c!(G)
    t_gen_cilk = @time_func generate_g_cilk!(G)
    t_gen_omp = @time_func generate_g_omp!(G)
    t_memset = @time_func fill_memset(G)
    t_fill = @time_func fill!(G, 0)
    t_fill_cilk = @time_func fill_cilk(G)
    t_fill_omp = @time_func fill_omp(G)

    return Dict("n" => n,
                "alloc" => t_alloc,
                "gen" => t_gen,
                "gen_c" => t_gen_c,
                "gen_cilk" => t_gen_cilk,
                "gen_omp" => t_gen_omp,
                "memset" => t_memset,
                "fill" => t_fill,
                "fill_cilk" => t_fill_cilk,
                "fill_omp" => t_fill_omp)
end

time_generates() =
    [time_generate(n) for n in [1, 2, 5, 10, 20, 50,
                                100, 200, 500, 1000, 2000, 5000, 10000]]

open(ARGS[1], "w") do io
    JSON.print(io, time_generates())
end
