using TemplateMatching
using Test
using StatsBase
using LinearAlgebra
using OffsetArrays
using CUDA

if CUDA.functional()
    @info "Running GPU tests" CUDA.version()
end

@testset verbose = true "TemplateMatching.jl" begin
    @testset "Cross-correlate" begin
        @test_throws ArgumentError crosscorrelate([1, 2, 3], Int[])
        @test_throws ArgumentError crosscorrelate(Int[], [1, 2, 3])
        @test_throws DimensionMismatch crosscorrelate([1, 2], [1, 2, 3])
        @test length(crosscorrelate([1, 2, 3], [1, 2])) == 2
        for usefft in (true, false), type = (Float32, Float64)
            let 
                x = @. 1e3 * ($rand(type, 100) - 0.5)
                y = @view x[45:55]
                y_mean, y_std = mean_and_std(y, corrected=false)
                z = @. (y - y_mean) / y_std
                @test all(crosscorrelate(x, z, type, normalize_template=false, usefft=usefft) .≈ crosscorrelate(x, y, type, usefft=usefft))
                if usefft && CUDA.functional()
                    x_d = CuArray(x)
                    z_d = CuArray(z)
                    @test all(crosscorrelate(x_d, z_d, type, normalize_template=false) .≈ crosscorrelate(x_d, z_d, type))
                end
            end
            let 
                x = @. 1e3 * ($rand(type, 100) - 0.5)
                y = @view x[45:55]
                cc, indx = findmax(crosscorrelate(x, y, type, usefft=usefft))
                @test cc ≈ 1.0 
                @test indx == 45
                if usefft && CUDA.functional()
                    x_d = CuArray(x)
                    y_d = CuArray(y)
                    cc, indx = findmax(crosscorrelate(x_d, y_d, type))
                    @test cc ≈ 1.0 
                    @test indx == 45
                end
            end
        end
    end

    @testset "Max filter" begin
        cases = [(Int[], 0, Int[]),
                 (Int[], 1, Int[]),
                 ([1, 2, 3], 0, [1, 2, 3]),
                 ([1, 2, 3], 1, [2, 3, 3]),
                 ([1, 2, 3], 2, [3, 3, 3])]
        for (u, l, v) in cases
            @test maxfilter(u, l) == v
        end
        if CUDA.functional()
            for (u, l, v) in cases
                @test maxfilter(CuArray(u), l) == CuArray(v)
            end
        end
    end

    @testset "Stack" begin
        @test_throws ArgumentError stack(Vector{Float64}[], Int[])
        @test_throws DimensionMismatch stack([[1, 2, 3], [4, 5, 6]], [1, 2, 3])
        @test stack([[I[i, j] ? 1 : 0 for i = 1:10] for j = 1:10], 1:10) == OffsetVector([1.0], 0:0)
        @test stack([[I[i, j] ? 10 : 0 for i = 1:10] for j = 1:10], zeros(Int, 10)) == ones(10)
    end

    @testset "Correlate template" begin
        @test_throws ArgumentError correlatetemplate(Vector{Float64}[], Vector{Float64}[], Int[], 0)
        @test_throws DimensionMismatch correlatetemplate([[1, 2, 3], [4, 5, 6], [7, 8, 9]], [[1, 2, 3], [4, 5, 6]], [1, 2], 0)
        for usefft in (true, false), type = (Float32, Float64)
            let 
                data = [rand(100) for _ = 1:10]
                template = [view(data[n], 45 + n:55 + n) for n = eachindex(data)]
                shifts = [45 + n for n = eachindex(data)]
                cc, indx = findmax(correlatetemplate(data, template, shifts, 0, type, usefft=usefft))
                @test isapprox(1.0, cc, atol=1e-5)
                @test indx == 0
                if usefft && CUDA.functional()
                    data_d = CuArray.(data)
                    template_d = CuArray.(template)
                    cc, indx = findmax(convert(OffsetVector{type, Vector{type}}, correlatetemplate(data_d, template_d, shifts, 0, type)))
                    @test isapprox(1.0, cc, atol=1e-5)
                    @test indx == 0
                end
            end
        end
    end

    @testset "Find Peaks" begin
        @test findpeaks([], 0, 0) == (Int[], [])
        let 
            x = rand(100)
            y = view(x, 40:60)
            cc = crosscorrelate(x, y)
            peaks, heights = findpeaks(cc, 0.99, 0)
            @test peaks == [40]
            @test heights ≈ [1.0]
        end
        let
            peaks, heights = findpeaks([x < pi ? 2 * abs(sin(x)) : abs(sin(x)) for x = 0:0.01pi:2pi], 0, 100)
            @test peaks == [51]
            @test heights ≈ [2.0]
        end
        let
            peaks, heights = findpeaks([x < pi ? 2 * abs(sin(x)) : abs(sin(x)) for x = 0:0.01pi:2pi], 1.0, 0)
            @test peaks == [51]
            @test heights ≈ [2.0]
        end
        let
            peaks, heights = findpeaks([x < pi ? abs(sin(x)) : 2 * abs(sin(x)) for x = 0:0.01pi:2pi], 1.0, 0)
            @test peaks == [151]
            @test heights ≈ [2.0]
        end
    end

    @testset "Relative magnitude" begin
        @test_throws ArgumentError relative_magnitude(Float64[], Float64[])
        @test_throws ArgumentError relative_magnitude([1, 2, 3], [1, 2])
        let num_trials = 1000, num_samples = 100, const_trace = fill(0.5, num_samples)
            @test isnan(relative_magnitude(ones(num_samples), zeros(num_samples)))
            @test isnan(relative_magnitude(zeros(num_samples), ones(num_samples)))
            @test relative_magnitude(const_trace, const_trace) == 0.0

            magnitude_test_func(factor) = mean(relative_magnitude(factor * rand(num_samples), rand(num_samples)) for _ = 1:num_trials)
            @test isapprox(magnitude_test_func(1.0), 0.0, atol=1e-2)
            @test isapprox(magnitude_test_func(10.0), 1.0, atol=1e-2)
        end
    end

    @testset "Find maximum in window" begin
        @test_throws ArgumentError findmax_window(Float64[], 0, 0)
        let x = rand(11)
            @test findmax_window(x, 6, 5) == findmax(x)
        end
    end

    @testset "Sub-sample shifting" begin
        let n_pts = 1024, tol = 1e-3
            x = sin.(range(0, 2pi, n_pts))
            for δ in [0.1, 0.5, 0.9] 
                x_reconstructed = TemplateMatching.subsampleshift(TemplateMatching.subsampleshift(x, δ), -δ)
                @test std(x .- x_reconstructed) < tol
            end
        end
    end
end
