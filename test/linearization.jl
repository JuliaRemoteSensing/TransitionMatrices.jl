@testitem "Linearization framework" begin
    problem = LinearizationProblem([2.0, 3.0, 1.311, 0.02, 2π];
                                   variables = (:a, :c, :mᵣ, :mᵢ, :λ)) do x
        (; shape = Spheroid(x[1], x[2], complex(x[3], x[4])), λ = x[5])
    end

    rebuilt = rebuild(problem)
    @test variables(problem) == (:a, :c, :mᵣ, :mᵢ, :λ)
    @test rebuilt.shape == Spheroid(2.0, 3.0, 1.311 + 0.02im)
    @test rebuilt.λ == 2π

    @test_throws ArgumentError LinearizationProblem([1.0, 2.0]; variables = (:a,)) do x
        x
    end

    scalar = LinearizationResult(10.0, [1.0, 2.0], (:x, :y))
    @test derivative(scalar, 1) == 1.0
    @test derivative(scalar, :y) == 2.0

    value = [1.0 2.0; 3.0 4.0]
    jacobian = reshape(collect(1.0:12.0), 2, 2, 3)
    matrix_result = LinearizationResult(value, jacobian, (:a, :c, :λ))
    @test derivative(matrix_result, :c) == view(jacobian, :, :, 2)

    support = supports_linearization(problem, EBCMLinearization())
    @test !Bool(support)
    @test !isempty(support.reason)
    @test_throws UnsupportedLinearization linearize_transition_matrix(problem, EBCMLinearization())
    @test_throws UnsupportedLinearization linearize_observable(scattering_cross_section,
                                                               problem,
                                                               IITMLinearization())
end

@testitem "Mie linearization matches ForwardDiff references" begin
    using ForwardDiff

    N = 5
    x₀ = [1.7, 1.311, 0.02, 2π]
    vars = (:x, :mᵣ, :mᵢ, :λ)

    problem = LinearizationProblem(x₀; variables = vars) do x
        (; x = x[1], m = complex(x[2], x[3]), λ = x[4], nₘₐₓ = N)
    end

    support = supports_linearization(problem, MieLinearization())
    @test Bool(support)

    result = linearize_transition_matrix(problem, MieLinearization())
    reference = MieTransitionMatrix{ComplexF64, N}(x₀[1], complex(x₀[2], x₀[3]))
    @test result.value.a ≈ reference.a
    @test result.value.b ≈ reference.b

    function coefficient_vector(x)
        T = eltype(x)
        a, b = TransitionMatrices.bhmie(T, x[1], complex(x[2], x[3]); nₘₐₓ = N)
        return vcat(real.(a), imag.(a), real.(b), imag.(b))
    end

    coefficient_jacobian = ForwardDiff.jacobian(coefficient_vector, x₀)
    for (j, var) in enumerate(vars)
        ∂T = derivative(result, var)
        ∂a = coefficient_jacobian[1:N, j] .+ coefficient_jacobian[(N + 1):(2N), j] .* im
        ∂b = coefficient_jacobian[(2N + 1):(3N), j] .+
             coefficient_jacobian[(3N + 1):(4N), j] .* im

        @test ∂T.a ≈ ∂a atol=1e-9 rtol=1e-9
        @test ∂T.b ≈ ∂b atol=1e-9 rtol=1e-9
    end

    function mie_scattering(x)
        T = eltype(x)
        𝐓 = MieTransitionMatrix{Complex{T}, N}(x[1], complex(x[2], x[3]))
        scattering_cross_section(𝐓, x[4])
    end

    function mie_extinction(x)
        T = eltype(x)
        𝐓 = MieTransitionMatrix{Complex{T}, N}(x[1], complex(x[2], x[3]))
        extinction_cross_section(𝐓, x[4])
    end

    function mie_absorption(x)
        T = eltype(x)
        𝐓 = MieTransitionMatrix{Complex{T}, N}(x[1], complex(x[2], x[3]))
        absorption_cross_section(𝐓, x[4])
    end

    function mie_albedo(x)
        T = eltype(x)
        𝐓 = MieTransitionMatrix{Complex{T}, N}(x[1], complex(x[2], x[3]))
        albedo(𝐓)
    end

    angles = (0.2, 0.3, 1.2, 0.5)
    function mie_amplitude_vector(x)
        T = eltype(x)
        𝐓 = MieTransitionMatrix{Complex{T}, N}(x[1], complex(x[2], x[3]))
        S = amplitude_matrix(𝐓, angles...; λ = x[4])
        return vcat(real.(vec(S)), imag.(vec(S)))
    end

    Csca = linearize_observable(scattering_cross_section, problem, MieLinearization())
    Cext = linearize_observable(extinction_cross_section, problem, MieLinearization())
    Cabs = linearize_observable(absorption_cross_section, problem, MieLinearization())
    ω = linearize_observable(albedo, problem, MieLinearization())
    S = linearize_observable(amplitude_matrix, problem, MieLinearization();
                             config = (; angles))

    @test Csca.value ≈ mie_scattering(x₀)
    @test Cext.value ≈ mie_extinction(x₀)
    @test Cabs.value ≈ mie_absorption(x₀)
    @test ω.value ≈ mie_albedo(x₀)
    @test S.value ≈ amplitude_matrix(reference, angles...; λ = x₀[4])
    @test Csca.jacobian ≈ ForwardDiff.gradient(mie_scattering, x₀) atol=1e-8 rtol=1e-8
    @test Cext.jacobian ≈ ForwardDiff.gradient(mie_extinction, x₀) atol=1e-8 rtol=1e-8
    @test Cabs.jacobian ≈ ForwardDiff.gradient(mie_absorption, x₀) atol=1e-8 rtol=1e-8
    @test ω.jacobian ≈ ForwardDiff.gradient(mie_albedo, x₀) atol=1e-8 rtol=1e-8

    amplitude_jacobian = ForwardDiff.jacobian(mie_amplitude_vector, x₀)
    for (j, var) in enumerate(vars)
        ∂S = derivative(S, var)
        ∂S_ref = reshape(amplitude_jacobian[1:4, j] .+
                         amplitude_jacobian[5:8, j] .* im, 2, 2)
        @test ∂S ≈ ∂S_ref atol=1e-8 rtol=1e-8
    end
end
