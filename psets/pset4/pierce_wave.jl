#!/usr/bin/julia -f

# We load the the data as UInt32 instead of UInt8, which is also enough
# to hold the bit map, since it generates more efficient code. Presumably
# because the CPU needs to do extra work to splat a <i8 x 8> to a 256bit ymm*
# register. Memory bandwidth doesn't seems to be an issue.
const map_data = readdlm("pierce_normalized.txt", ' ', UInt32)
const orig_map = readdlm("pierce.txt", ' ', UInt32)

@inline color_map_r(v::Float32) =
    ifelse(v < 0.5f0, ifelse(v < 0.25f0, 1f0, (0.5f0 - v) * 4), 0f0)
@inline color_map_g(v::Float32) =
    ifelse(v < 0.75f0, ifelse(v < 0.25f0, v * 4, 1f0), (1f0 - v) * 4)
@inline color_map_b(v::Float32) =
    ifelse(v < 0.5f0, 0f0, ifelse(v < 0.75f0, (v - 0.5f0) * 4, 1f0))

immutable ZeroInitializer
end

function Base.call(::ZeroInitializer, map_data, prev_buff, cur_buff)
    fill!(prev_buff, 0)
    fill!(cur_buff, 0)
end

immutable SinDrive
    p_0::Float32
    ω::Float32
    xrng::UnitRange{Int}
    yrng::UnitRange{Int}
end

function Base.call(d::SinDrive, buff, t)
    v = d.p_0 * sin(d.ω * t)
    @inbounds for j in d.yrng
        for i in d.xrng
            buff[i, j] = v
        end
    end
    nothing
end

immutable DummyCollector
end

Base.call(::DummyCollector, buff, idx) = nothing

immutable FullCollector
    result::Array{Float32,4}
    skip::Int
    FullCollector(map_data, nsteps, skip=1) =
        new(Array{Float32}(size(map_data)..., 3,
                           nsteps ÷ skip), skip)
end

function Base.call(c::FullCollector, buff, idx)
    (idx - 1) % c.skip != 0 && return
    idx = (idx - 1) ÷ c.skip + 1
    res = c.result
    @assert idx <= size(res, 4)
    vmax = max(-minimum(buff), maximum(buff))
    @inbounds for j in 1:size(res, 2)
        for i in 1:size(res, 1)
            px_t = map_data[i, j]
            if px_t != 0
                v = buff[i, j] / vmax
                v = flipsign(v, abs(v)^0.05f0)
                v = (v + 1) / 2
                res[i, j, 1, idx] = color_map_r(v)
                res[i, j, 2, idx] = color_map_g(v)
                res[i, j, 3, idx] = color_map_b(v)
            elseif orig_map[i, j] != 0
                res[i, j, 1, idx] = 0
                res[i, j, 2, idx] = 0
                res[i, j, 3, idx] = 0
            else
                res[i, j, 1, idx] = 1
                res[i, j, 2, idx] = 1
                res[i, j, 3, idx] = 1
            end
        end
    end
    nothing
end

function propagate_wave(map_data, initializer, drive_override, c, _Δt, h,
                        nsteps, result_collector)
    # Initialize some useful constants
    Δt = Float32(_Δt)
    # Effective c is c / h * Δt
    c_eff = Float32(c / h * Δt)
    c2 = c_eff^2
    #   Pⁿ⁺¹ⱼₖ
    # = 2Pⁿⱼₖ - Pⁿ⁻¹ⱼₖ + c2 * (Pⁿⱼ₊₁ₖ + Pⁿⱼₖ₊₁ + Pⁿⱼ₋₁ₖ + Pⁿⱼₖ₋₁ - 4Pⁿⱼₖ)
    # = (2 - 4 * c2) * Pⁿⱼₖ - Pⁿ⁻¹ⱼₖ + c2 * (Pⁿⱼ₊₁ₖ + Pⁿⱼₖ₊₁ + Pⁿⱼ₋₁ₖ + Pⁿⱼₖ₋₁)
    coeff_jk = 2 - 4 * c2
    # = muladd(c2, (Pⁿⱼ₊₁ₖ + Pⁿⱼₖ₊₁) + (Pⁿⱼ₋₁ₖ + Pⁿⱼₖ₋₁),
    #          muladd(coeff_jk, Pⁿⱼₖ, -Pⁿ⁻¹ⱼₖ))

    # Initialize our buffers
    xlen, ylen = size(map_data)
    prev_buff = Matrix{Float32}(xlen, ylen)
    cur_buff = Matrix{Float32}(xlen, ylen)
    # Fill them with the initializer supplied
    initializer(map_data, prev_buff, cur_buff)

    drive_override(prev_buff, Δt)
    result_collector(prev_buff, 1)
    drive_override(prev_buff, 2Δt)
    result_collector(cur_buff, 2)

    # Now we can determine indexing range we actualy want to loop over in each
    # column. In principle, we should divide the range into sub ranges when
    # there's a big gap in between. However, in the direction we are looping
    # over in the inner loop, the largest gap is ~ 8 pixels, which is basically
    # the vector size for f32 with AVX2 instructions and smaller than the size
    # of the cache line (cache line might not matter too much though since the
    # whole buffer fits in L2 cache). Therefore, we simply figure out the lower
    # and upper bound for each index. (A more general purpose function should
    # record the column number, lower bound, upper bound triplet instead).
    min_yidx = 1
    max_yidx = 0
    idx_ranges = Vector{UnitRange{Int}}(ylen)
    # The 2:(len - 1) loop range is a really minor optimization since we
    # already assume the edges are not interesting.
    idx_ranges[1] = 1:0
    @inbounds for j in 2:(ylen - 1)
        min_xidx = 1
        max_xidx = 0
        for i in 2:(xlen - 1)
            if map_data[i, j] != 0
                max_xidx == 0 && (min_xidx = i)
                max_xidx = i
            end
        end
        if max_xidx != 0
            max_yidx == 0 && (min_yidx = j)
            max_yidx = j
        end
        idx_ranges[j] = min_xidx:max_xidx
    end
    idx_ranges[end] = 1:0

    @inbounds for t in 3:nsteps
        # prev_buff serve both as the previous buffer and the next buffer
        # since we only read the current element form it once
        @fastmath for j in min_yidx:max_yidx
            @simd for i in idx_ranges[j]
                # Load memories
                v = cur_buff[i, j]
                v_tm1 = -prev_buff[i, j]
                v_ym1 = cur_buff[i, j - 1]
                v_xm1 = cur_buff[i - 1, j]
                v_xp1 = cur_buff[i + 1, j]
                v_yp1 = cur_buff[i, j + 1]
                px_t = map_data[i, j]

                # Replace wall value with self value
                v_ym1 = ifelse(px_t & 0x2 != 0, v_ym1, v)
                v_xm1 = ifelse(px_t & 0x4 != 0, v_xm1, v)
                v_xp1 = ifelse(px_t & 0x8 != 0, v_xp1, v)
                v_yp1 = ifelse(px_t & 0x10 != 0, v_yp1, v)

                m1 = muladd(coeff_jk, v, v_tm1)
                vsum = (v_ym1 + v_xm1) + (v_yp1 + v_xp1)
                prev_buff[i, j] = muladd(c2, vsum, m1)
            end
        end

        prev_buff, cur_buff = cur_buff, prev_buff
        drive_override(cur_buff, t * Δt)
        result_collector(cur_buff, t)
    end
    nothing
end

immutable FramesCollector{N}
    frames::NTuple{N,Array{Float32,3}}
    frame_nums::NTuple{N,Int}
    @generated function FramesCollector(map_data, frame_nums::NTuple{N,Int})
        quote
            map_size = size(map_data)
            frames = ($([:(Array{Float32}(map_size..., 3)) for i in 1:N]...),)
            $(Expr(:new, FramesCollector{N}, :frames, :frame_nums))
        end
    end
end

FramesCollector{N}(map_data, frame_nums::NTuple{N,Int}) =
    FramesCollector{N}(map_data, frame_nums)

function Base.call{N}(c::FramesCollector{N}, buff, idx)
    frame_nums = c.frame_nums
    @inbounds for k in 1:N
        # O(N)... whatever...
        if idx == frame_nums[k]
            vmax = max(-minimum(buff), maximum(buff))
            frame = c.frames[k]
            for j in 1:size(frame, 2)
                for i in 1:size(frame, 1)
                    px_t = map_data[i, j]
                    if px_t != 0
                        v = (buff[i, j] / vmax + 1) / 2
                        frame[i, j, 1] = color_map_r(v)
                        frame[i, j, 2] = color_map_g(v)
                        frame[i, j, 3] = color_map_b(v)
                    elseif orig_map[i, j] != 0
                        frame[i, j, 1] = 0
                        frame[i, j, 2] = 0
                        frame[i, j, 3] = 0
                    else
                        frame[i, j, 1] = 1
                        frame[i, j, 2] = 1
                        frame[i, j, 3] = 1
                    end
                end
            end
        end
    end
    nothing
end

function main(nstep, Δt)
    # collector = DummyCollector()
    collector = FullCollector(map_data, nstep, 5)
    # collector = FramesCollector(map_data, (Int(0.015 ÷ Δt) + 1,
    #                                        Int(0.105 ÷ Δt) + 1,
    #                                        Int(0.505 ÷ Δt) + 1,
    #                                        Int(1.005 ÷ Δt) + 1))
    @time propagate_wave(map_data, ZeroInitializer(),
                         SinDrive(10f0, 100f0π, 58:61, 16:19),
                         3.43f4, Δt, 36.6f0, nstep, collector)
    collector
end
# const collector = main(20_000, 1f-4)

# using PyPlot
# figure()
# imshow(collector.frames[1], interpolation="none")
# title("0.015s")
# savefig("pierce_015")
# figure()
# imshow(collector.frames[2], interpolation="none")
# title("0.105s")
# savefig("pierce_105")
# figure()
# imshow(collector.frames[3], interpolation="none")
# title("0.505s")
# savefig("pierce_505")
# figure()
# imshow(collector.frames[4], interpolation="none")
# title("1.005s")
# savefig("pierce_1005")
# show()

# The following rather ugly code is for initial benchmarking and making videos

# const collector = main(20_000, 1f-4)

# using PyCall
# using PyPlot
# const animation = pyimport("matplotlib.animation")

# const fig = figure(figsize=(100, 200))
# const ax = gca()
# ax[:set_xticks]([])
# ax[:set_yticks]([])

# const im = imshow(collector.result[:, :, :, 1], interpolation="none")

# # initialization function: plot the background of each frame
# plot_init() = plot_animate(0)

# # animation function.  This is called sequentially
# function plot_animate(i)
#     im[:set_data](collector.result[:, :, :, i + 1])
#     (im,)
# end

# @time anim = animation[:FuncAnimation](fig, plot_animate, init_func=plot_init,
#                                        frames=size(collector.result, 4),
#                                        interval=16)
# # Somehow the saver only saves the part of the animation shown on the screen...
# show()
# FFwriter = animation[:FFMpegWriter](fps=60)
# @time anim[:save]("pierce.mp4", writer=FFwriter, fps=60,
#                   extra_args=["-vcodec", "libx264"])
