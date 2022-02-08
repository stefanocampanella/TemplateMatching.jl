using TemplateMatching
using Test
using StatsBase
using LinearAlgebra
using OffsetArrays



@testset verbose = true "TemplateMatching.jl" begin
    @testset "Cross-correlate" begin
        @test_throws ArgumentError crosscorrelate([1, 2, 3], Int[])
        @test_throws ArgumentError crosscorrelate(Int[], [1, 2, 3])
        @test_throws DimensionMismatch crosscorrelate([1, 2], [1, 2, 3])
        @test length(crosscorrelate([1, 2, 3], [1, 2])) == 2
        let 
            x = @. 1e3 * ($rand(100) - 0.5)
            y = @view x[45:55]
            y_mean, y_std = mean_and_std(y, corrected=false)
            z = @. (y - y_mean) / y_std
            @test all(crosscorrelate(x, z, normalize_template=false) .≈ crosscorrelate(x, y) )
        end
        let 
            x = @. 1e3 * ($rand(100) - 0.5)
            y = @view x[45:55]
            cc, indx = findmax(crosscorrelate(x, y))
            @test cc ≈ 1.0 
            @test indx == 45
        end
    end

    @testset "Max filter" begin
        @test maxfilter([], 0) == []
        @test maxfilter([], 1) == []
        @test maxfilter([], 2) == []
        @test maxfilter([1, 2, 3], 0) == [1, 2, 3]
        @test maxfilter([1, 2, 3], 1) == [2, 3, 3]
        @test maxfilter([1, 2, 3], 2) == [3, 3, 3]
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
        let 
            data = [rand(100) for _ = 1:10]
            template = [view(data[n], 45 + n:55 + n) for n = eachindex(data)]
            shifts = [45 + n for n = eachindex(data)]
            cc, indx = findmax(correlatetemplate(data, template, shifts, 0))
            @test cc ≈ 1.0 
            @test indx == 0
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
end
