using TemplateMatching
using Test
using StatsBase
using LinearAlgebra


@testset verbose = true "TemplateMatching.jl" begin
    @testset "Cross-correlate" begin
        @test_throws ArgumentError crosscorrelate([1, 2, 3], Int[])
        @test_throws ArgumentError crosscorrelate(Int[], [1, 2, 3])
        @test_throws DimensionMismatch crosscorrelate([1, 2], [1, 2, 3])
        @test length(crosscorrelate([1, 2, 3], [1, 2])) == 2
        @test let 
            x = @. 1e3 * ($rand(100) - 0.5)
            y = @view x[45:55]
            y_mean, y_std = mean_and_std(y, corrected=false)
            z = @. (y - y_mean) / y_std
            all(crosscorrelate(x, z, normalize_template=false) .≈ crosscorrelate(x, y) )
        end
        @test let 
            x = @. 1e3 * ($rand(100) - 0.5)
            y = @view x[45:55]
            cc, indx = findmax(crosscorrelate(x, y))
            cc ≈ 1.0 && indx == 45
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
        @test stack([[I[i, j] ? 1 : 0 for i = 1:10] for j = 1:10], 1:10) == [1.0]
        @test stack([[I[i, j] ? 10 : 0 for i = 1:10] for j = 1:10], zeros(Int, 10)) == ones(10)

    end
end
