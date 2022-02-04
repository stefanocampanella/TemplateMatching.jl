using TemplateMatching
using Test
using StatsBase

@testset "TemplateMatching.jl" begin
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
