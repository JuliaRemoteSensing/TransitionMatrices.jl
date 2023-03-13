@testitem "Auto differentiation" begin
    using ForwardDiff, FiniteDifferences

    function f(x)
        s = Spheroid(x[1], x[2], complex(x[3], x[4]))
        T₀ = TransitionMatrices.transition_matrix(s, x[5], 10, 100)
        Csca = calc_Csca(T₀)
    end

    gradient = ForwardDiff.gradient(f, [2.0, 3.0, 1.311, 0.02, 2π])

    @test gradient[1]≈central_fdm(5, 1)(x -> f([x, 3.0, 1.311, 0.02, 2π]), 2.0) atol=1e-6
    @test gradient[2]≈central_fdm(5, 1)(x -> f([2.0, x, 1.311, 0.02, 2π]), 3.0) atol=1e-6
    @test gradient[3]≈central_fdm(5, 1)(x -> f([2.0, 3.0, x, 0.02, 2π]), 1.311) atol=1e-6
    @test gradient[4]≈central_fdm(5, 1)(x -> f([2.0, 3.0, 1.311, x, 2π]), 0.02) atol=1e-6
    @test gradient[5]≈central_fdm(5, 1)(x -> f([2.0, 3.0, 1.311, 0.02, x]), 2π) atol=1e-6
end
