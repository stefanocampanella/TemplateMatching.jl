using TemplateMatching
using Test

@testset "TemplateMatching.jl" begin
    @test_throws ErrorException crosscorrelate([1, 2, 3], Int[])
    @test_throws ErrorException crosscorrelate(Int[], [1, 2, 3])
    @test length(crosscorrelate([1, 2, 3], [1, 2])) == 2
    @test let 
        x = @. 100.0 * ($rand(100) - 0.5)
        y = @view x[45:55]
        cc, indx = findmax(crosscorrelate(x, y))
        cc ≈ 1.0 && indx == 45
    end
end
