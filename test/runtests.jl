using TransitionMatrices
using Test
using TestItemRunner

@testset "TransitionMatrices.jl" begin
    @testset "Unit tests" begin @run_package_tests() end

    @testset "Integration tests" begin include("autodiff.jl") end
end
