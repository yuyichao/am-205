#!/usr/bin/julia -f

immutable HeatPropagator{T}
    α::Float32
    r_1::Float32
    nr::Int
    Δr::Float32
    Δt::Float32
    solver::T
    HeatPropagator(α, r_1, nr, Δr, Δt, solver) =
        new(α, r_1, nr, Δr, Δt, solver)
end

function HeatPropagator(α, r_1, nr, Δr, Δt, boundary)
    solver_orig = zeros(Float32, (nr, nr))
    @inbounds @fastmath for j in 2:(nr - 1)
        r = r_1 + (j - 1) * Δr
        solver_orig[j, j + 1] = -α * Δt / Δr * (1 / Δr + 1 / 2 / r)
        solver_orig[j, j] = 2 * α * Δt / Δr^2 + 1
        solver_orig[j, j - 1] = -α * Δt / Δr * (1 / Δr - 1 / 2 / r)
    end
    boundary(solver_orig, α, r_1, nr, Δr, Δt)
    solver = factorize(solver_orig)
    HeatPropagator{typeof(solver)}(α, r_1, nr, Δr, Δt, solver)
end

immutable OnesInitializer
end

function Base.call(::OnesInitializer, prop, buff)
    fill!(buff, 1)
end

immutable BoundaryP3
end

function Base.call(d::BoundaryP3, solver_orig, α, r_1, nr, Δr, Δt)
    solver_orig[1, 1] = 2 * α * Δt / Δr^2 + 1
    solver_orig[1, 2] = -2 * α * Δt / Δr^2
    r_n = r_1 + (nr - 1) * Δr
    solver_orig[nr, nr - 1] = -2 * α * Δt / Δr^2
    solver_orig[nr, nr] = 2 * α * Δt / Δr^2 + 1 + 2 * Δt / Δr + Δt / r_n
end

immutable DummyCollector
end

Base.call(::DummyCollector, buff, prop, idx) = nothing

function propagate_heat(initializer, propagator, nsteps, result_collector)
    buff = Vector{Float32}(propagator.nr)
    # Fill them with the initializer supplied
    initializer(propagator, buff)
    result_collector(buff, propagator, 1)

    solver = propagator.solver
    @inbounds @fastmath for i in 2:nsteps
        A_ldiv_B!(solver, buff)
        result_collector(buff, propagator, i)
    end
    nothing
end

immutable FramesCollector{N}
    frames::NTuple{N,Vector{Float32}}
    frame_nums::NTuple{N,Int}
    @generated function FramesCollector(prop, frame_nums::NTuple{N,Int})
        quote
            frames = ($([:(Vector{Float32}(prop.nr)) for i in 1:N]...),)
            $(Expr(:new, FramesCollector{N}, :frames, :frame_nums))
        end
    end
end

FramesCollector{N}(prop, frame_nums::NTuple{N,Int}) =
    FramesCollector{N}(prop, frame_nums)

function Base.call{N}(c::FramesCollector{N}, buff, prop, idx)
    frame_nums = c.frame_nums
    @inbounds for i in 1:N
        # O(N)... whatever...
        if idx == frame_nums[i]
            frame = c.frames[i]
            @simd for j in 1:prop.nr
                frame[j] = buff[j]
            end
            return
        end
    end
    nothing
end
