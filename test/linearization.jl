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
    @test derivative(matrix_result, :c) == view(jacobian,:,:,2)
    @test_throws BoundsError derivative(matrix_result, 0)
    @test_throws BoundsError derivative(matrix_result, 4)
    @test_throws ArgumentError derivative(matrix_result, :unknown)

    support = supports_linearization(problem, EBCMLinearization())
    @test !Bool(support)
    @test !isempty(support.reason)
    err = try
        linearize_transition_matrix(problem, EBCMLinearization())
    catch caught
        caught
    end
    @test err isa UnsupportedLinearization
    @test occursin("Unsupported linearization", sprint(showerror, err))
    @test_throws UnsupportedLinearization linearize_observable(scattering_cross_section,
        problem,
        IITMLinearization())
    @test_throws UnsupportedLinearization linearize_observable(:scattering_cross_section,
        problem,
        EBCMLinearization())
    @test_throws ArgumentError LinearizationProblem([1.0]; variables = ("x",)) do x
        x
    end
end

@testitem "IITM linearization support contracts" begin
    base_problem = LinearizationProblem([1.311]; variables = (:mᵣ,)) do x
        (; shape = Spheroid(0.9, 1.2, complex(x[1], 0.02)), λ = 2π)
    end

    missing_config = supports_linearization(base_problem, IITMLinearization())
    @test !Bool(missing_config)
    @test occursin("requires shape", missing_config.reason)

    config = (; nₘₐₓ = 2, Nr = 3, Nϑ = 8)
    @test Bool(supports_linearization(base_problem, IITMLinearization(); config))
    @test Bool(supports_linearization(base_problem, IITMLinearization(:auto); config))
    @test Bool(supports_linearization(base_problem, IITMLinearization(:axisymmetric);
        config))

    unsupported_output = supports_linearization(base_problem,
        IITMLinearization(:axisymmetric);
        output = :observable, config)
    @test !Bool(unsupported_output)
    @test occursin("only supports transition matrices", unsupported_output.reason)

    bad_variant = supports_linearization(base_problem, IITMLinearization(:bad); config)
    @test !Bool(bad_variant)
    @test occursin("variant must be", bad_variant.reason)

    geometry_problem = LinearizationProblem([0.9]; variables = (:a,)) do x
        (; shape = Spheroid(x[1], 1.2, 1.311 + 0.02im), λ = 2π)
    end
    geometry_support = supports_linearization(geometry_problem,
        IITMLinearization(:axisymmetric);
        config)
    @test !Bool(geometry_support)
    @test occursin("fixed-geometry", geometry_support.reason)

    duplicate_problem = LinearizationProblem([1.311, 1.4];
        variables = (:mᵣ, :mᵣ)) do x
        (; shape = Spheroid(0.9, 1.2, complex(x[1], 0.02)), λ = 2π)
    end
    duplicate_support = supports_linearization(duplicate_problem,
        IITMLinearization(:axisymmetric);
        config)
    @test !Bool(duplicate_support)
    @test occursin("unique canonical variables", duplicate_support.reason)

    prism_problem = LinearizationProblem([1.311]; variables = (:mᵣ,)) do x
        (; shape = Prism(5, 0.8, 1.1, complex(x[1], 0.02)), λ = 2π)
    end
    axisymmetric_support = supports_linearization(prism_problem,
        IITMLinearization(:axisymmetric);
        config)
    @test !Bool(axisymmetric_support)
    @test occursin("requires an axisymmetric shape", axisymmetric_support.reason)

    nfold_missing_phi = supports_linearization(prism_problem,
        IITMLinearization(:nfold); config)
    @test !Bool(nfold_missing_phi)
    @test occursin("requires Nφ", nfold_missing_phi.reason)

    arbitrary_missing_phi = supports_linearization(prism_problem,
        IITMLinearization(:arbitrary);
        config)
    @test !Bool(arbitrary_missing_phi)
    @test occursin("requires Nφ", arbitrary_missing_phi.reason)
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

        @test ∂T.a≈∂a atol=1e-9 rtol=1e-9
        @test ∂T.b≈∂b atol=1e-9 rtol=1e-9
    end

    subset_x₀ = [x₀[2], x₀[1]]
    subset_vars = (:mᵣ, :x)
    subset_problem = LinearizationProblem(subset_x₀; variables = subset_vars) do x
        (; x = x[2], m = complex(x[1], x₀[3]), λ = x₀[4], nₘₐₓ = N)
    end

    function subset_coefficient_vector(x)
        T = eltype(x)
        a, b = TransitionMatrices.bhmie(T, x[2], complex(x[1], x₀[3]); nₘₐₓ = N)
        return vcat(real.(a), imag.(a), real.(b), imag.(b))
    end

    @test Bool(supports_linearization(subset_problem, MieLinearization()))
    subset_result = linearize_transition_matrix(subset_problem, MieLinearization())
    subset_coefficient_jacobian = ForwardDiff.jacobian(subset_coefficient_vector, subset_x₀)

    for (j, var) in enumerate(subset_vars)
        ∂T = derivative(subset_result, var)
        ∂a = subset_coefficient_jacobian[1:N, j] .+
             subset_coefficient_jacobian[(N + 1):(2N), j] .* im
        ∂b = subset_coefficient_jacobian[(2N + 1):(3N), j] .+
             subset_coefficient_jacobian[(3N + 1):(4N), j] .* im

        @test ∂T.a≈∂a atol=1e-9 rtol=1e-9
        @test ∂T.b≈∂b atol=1e-9 rtol=1e-9
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
    @test Csca.jacobian≈ForwardDiff.gradient(mie_scattering, x₀) atol=1e-8 rtol=1e-8
    @test Cext.jacobian≈ForwardDiff.gradient(mie_extinction, x₀) atol=1e-8 rtol=1e-8
    @test Cabs.jacobian≈ForwardDiff.gradient(mie_absorption, x₀) atol=1e-8 rtol=1e-8
    @test ω.jacobian≈ForwardDiff.gradient(mie_albedo, x₀) atol=1e-8 rtol=1e-8

    amplitude_jacobian = ForwardDiff.jacobian(mie_amplitude_vector, x₀)
    for (j, var) in enumerate(vars)
        ∂S = derivative(S, var)
        ∂S_ref = reshape(
            amplitude_jacobian[1:4, j] .+
            amplitude_jacobian[5:8, j] .* im, 2, 2)
        @test ∂S≈∂S_ref atol=1e-8 rtol=1e-8
    end
end

@testitem "EBCM block linearization follows the P-U matrix identity" begin
    using TransitionMatrices: 𝐓_from_𝐏_and_𝐔, ∂𝐓_from_𝐏_and_𝐔

    𝐏 = ComplexF64[1.2 + 0.1im 0.05 - 0.03im
                   -0.04 + 0.02im 0.9 + 0.2im]
    𝐔 = ComplexF64[0.4 - 0.2im 0.08 + 0.04im
                   0.03 - 0.06im 0.5 + 0.1im]
    ∂𝐏 = ComplexF64[0.2 - 0.1im -0.04 + 0.03im
                    0.01 + 0.05im -0.08 + 0.02im]
    ∂𝐔 = ComplexF64[-0.1 + 0.03im 0.07 - 0.01im
                    -0.02 + 0.04im 0.06 + 0.05im]

    ∂𝐓 = ∂𝐓_from_𝐏_and_𝐔(𝐏, 𝐔, ∂𝐏, ∂𝐔)
    ϵ = 1e-6
    ∂𝐓_fd = (𝐓_from_𝐏_and_𝐔(𝐏 .+ ϵ .* ∂𝐏, 𝐔 .+ ϵ .* ∂𝐔) -
             𝐓_from_𝐏_and_𝐔(𝐏 .- ϵ .* ∂𝐏, 𝐔 .- ϵ .* ∂𝐔)) ./ (2ϵ)

    @test ∂𝐓≈∂𝐓_fd atol=1e-9 rtol=1e-8
end

@testitem "EBCM block assembly exposes P-U matrices" begin
    using TransitionMatrices: Spheroid, 𝐓_from_𝐏_and_𝐔, ebcm_matrices_m,
                              ebcm_matrices_m₀, transition_matrix_m,
                              transition_matrix_m₀

    s = Spheroid{Float64, ComplexF64}(1.0, 1.2, 1.311 + 0.02im)
    λ = 2π
    nₘₐₓ = 4
    Ng = 40

    𝐏₀, 𝐔₀ = ebcm_matrices_m₀(s, λ, nₘₐₓ, Ng)
    @test 𝐓_from_𝐏_and_𝐔(𝐏₀, 𝐔₀) ≈ transition_matrix_m₀(s, λ, nₘₐₓ, Ng)

    m = 2
    𝐏ₘ, 𝐔ₘ = ebcm_matrices_m(m, s, λ, nₘₐₓ, Ng)
    @test 𝐓_from_𝐏_and_𝐔(𝐏ₘ, 𝐔ₘ) ≈ transition_matrix_m(m, s, λ, nₘₐₓ, Ng)
end

@testitem "EBCM block matrices assemble axisymmetric transition linearization" begin
    using TransitionMatrices: ebcm_transition_matrix_from_matrices,
                              ∂ebcm_transition_matrix_from_matrices,
                              𝐓_from_𝐏_and_𝐔, ∂𝐓_from_𝐏_and_𝐔

    function test_matrix(n, offset; scale = 1.0)
        [scale * ((i == j ? 1.5 + offset : 0.03 * (i + j + offset)) +
          0.02im * (i - j + offset)) for i in 1:n, j in 1:n]
    end

    𝐏s = [test_matrix(4, 1), test_matrix(4, 2), test_matrix(2, 3)]
    𝐔s = [test_matrix(4, 4; scale = 0.2),
        test_matrix(4, 5; scale = 0.2),
        test_matrix(2, 6; scale = 0.2)]
    ∂𝐏s = [test_matrix(4, 7; scale = 0.1),
        test_matrix(4, 8; scale = 0.1),
        test_matrix(2, 9; scale = 0.1)]
    ∂𝐔s = [test_matrix(4, 10; scale = 0.1),
        test_matrix(4, 11; scale = 0.1),
        test_matrix(2, 12; scale = 0.1)]

    𝐓 = ebcm_transition_matrix_from_matrices(𝐏s, 𝐔s)
    ∂𝐓 = ∂ebcm_transition_matrix_from_matrices(𝐏s, 𝐔s, ∂𝐏s, ∂𝐔s)

    for m in 0:2
        @test 𝐓.𝐓[m + 1] ≈ 𝐓_from_𝐏_and_𝐔(𝐏s[m + 1], 𝐔s[m + 1])
        @test ∂𝐓.𝐓[m + 1] ≈ ∂𝐓_from_𝐏_and_𝐔(𝐏s[m + 1], 𝐔s[m + 1],
            ∂𝐏s[m + 1], ∂𝐔s[m + 1])
    end
end

@testitem "IITM transition recurrence solve form matches finite differences" begin
    using TransitionMatrices: _iitm_update_transition_solve_block

    function test_matrix(n, offset; scale = 1.0)
        [scale * ((i == j ? 1.3 + offset : 0.02 * (i + j + offset)) +
          0.03im * (i - j + offset)) for i in 1:n, j in 1:n]
    end

    𝐓 = test_matrix(4, 1; scale = 0.1)
    𝐐ⱼⱼ = test_matrix(4, 2; scale = 0.08)
    𝐐ⱼₕ = test_matrix(4, 3; scale = 0.03)
    𝐐ₕⱼ = test_matrix(4, 4; scale = 0.03)
    𝐐ₕₕ = test_matrix(4, 5; scale = 0.04)
    ∂𝐓 = test_matrix(4, 6; scale = 0.01)
    ∂𝐐ⱼⱼ = test_matrix(4, 7; scale = 0.01)
    ∂𝐐ⱼₕ = test_matrix(4, 8; scale = 0.01)
    ∂𝐐ₕⱼ = test_matrix(4, 9; scale = 0.01)
    ∂𝐐ₕₕ = test_matrix(4, 10; scale = 0.01)

    𝐓next,
    ∂𝐓next = _iitm_update_transition_solve_block(𝐓, ∂𝐓, 𝐐ⱼⱼ,
        𝐐ⱼₕ, 𝐐ₕⱼ, 𝐐ₕₕ,
        ∂𝐐ⱼⱼ, ∂𝐐ⱼₕ,
        ∂𝐐ₕⱼ, ∂𝐐ₕₕ)

    update_value(ϵ) = _iitm_update_transition_solve_block(𝐓 .+ ϵ .* ∂𝐓, zero(∂𝐓),
        𝐐ⱼⱼ .+ ϵ .* ∂𝐐ⱼⱼ,
        𝐐ⱼₕ .+ ϵ .* ∂𝐐ⱼₕ,
        𝐐ₕⱼ .+ ϵ .* ∂𝐐ₕⱼ,
        𝐐ₕₕ .+ ϵ .* ∂𝐐ₕₕ,
        zero(∂𝐐ⱼⱼ), zero(∂𝐐ⱼₕ),
        zero(∂𝐐ₕⱼ), zero(∂𝐐ₕₕ))[1]

    ϵ = 1e-6
    ∂𝐓_fd = (update_value(ϵ) - update_value(-ϵ)) ./ (2ϵ)
    @test 𝐓next ≈ update_value(0)
    @test ∂𝐓next≈∂𝐓_fd atol=1e-8 rtol=1e-7
end

@testitem "EBCM fixed spheroid linearization matches ForwardDiff references" begin
    using ForwardDiff

    nₘₐₓ = 3
    Ng = 24
    x₀ = [1.0, 1.2, 1.311, 0.02, 2π]
    vars = (:a, :c, :mᵣ, :mᵢ, :λ)
    config = (; nₘₐₓ, Ng)

    problem = LinearizationProblem(x₀; variables = vars) do x
        (; shape = Spheroid(x[1], x[2], complex(x[3], x[4])), λ = x[5])
    end

    function transition_vector_from_parameters(x)
        s = Spheroid(x[1], x[2], complex(x[3], x[4]))
        𝐓 = TransitionMatrices.transition_matrix(s, x[5], nₘₐₓ, Ng)
        chunks = [vcat(real.(vec(block)), imag.(vec(block))) for block in 𝐓.𝐓]
        return reduce(vcat, chunks)
    end

    function transition_vector_from_matrix(𝐓)
        chunks = [vcat(real.(vec(block)), imag.(vec(block))) for block in 𝐓.𝐓]
        return reduce(vcat, chunks)
    end

    support = supports_linearization(problem, EBCMLinearization(); config)
    @test Bool(support)

    result = linearize_transition_matrix(problem, EBCMLinearization(); config)
    @test result.metadata.backend == :ebcm_analytic
    @test !hasproperty(result.metadata, :reference)
    reference = TransitionMatrices.transition_matrix(
        Spheroid(x₀[1], x₀[2],
            complex(x₀[3], x₀[4])),
        x₀[5], nₘₐₓ, Ng)
    reference_jacobian = ForwardDiff.jacobian(transition_vector_from_parameters, x₀)

    @test transition_vector_from_matrix(result.value) ≈
          transition_vector_from_matrix(reference)

    for (j, var) in enumerate(vars)
        @test transition_vector_from_matrix(derivative(result, var))≈reference_jacobian[:, j] atol=1e-7 rtol=1e-7
    end

    subset_x₀ = [x₀[5], x₀[1]]
    subset_vars = (:λ, :a)
    subset_problem = LinearizationProblem(subset_x₀; variables = subset_vars) do x
        c = zero(x[1]) + zero(x[2]) + x₀[2]
        (; shape = Spheroid(x[2], c, complex(x₀[3], x₀[4])),
            λ = x[1])
    end

    function transition_vector_from_subset(x)
        base = zero(x[1]) + zero(x[2])
        s = Spheroid(base + x[2], base + x₀[2],
            complex(base + x₀[3], base + x₀[4]))
        𝐓 = TransitionMatrices.transition_matrix(s, x[1], nₘₐₓ, Ng)
        return transition_vector_from_matrix(𝐓)
    end

    @test Bool(supports_linearization(subset_problem, EBCMLinearization(); config))
    subset_result = linearize_transition_matrix(subset_problem, EBCMLinearization(); config)
    subset_jacobian = ForwardDiff.jacobian(transition_vector_from_subset, subset_x₀)

    for (j, var) in enumerate(subset_vars)
        @test transition_vector_from_matrix(derivative(subset_result, var))≈subset_jacobian[:, j] atol=1e-7 rtol=1e-7
    end

    real_m_problem = LinearizationProblem([x₀[5]]; variables = (:λ,)) do x
        (; shape = Spheroid(x₀[1], x₀[2], x₀[3]), λ = x[1])
    end
    @test Bool(supports_linearization(real_m_problem, EBCMLinearization(); config))
    real_m_result = linearize_transition_matrix(real_m_problem, EBCMLinearization(); config)
    @test real_m_result.metadata.backend == :ebcm_analytic
end

@testitem "EBCM fixed Chebyshev linearization matches ForwardDiff references" begin
    using ForwardDiff

    nₘₐₓ = 3
    Ng = 28
    degree = 4
    x₀ = [1.0, 0.08, 1.311, 0.02, 2π]
    vars = (:r₀, :ε, :mᵣ, :mᵢ, :λ)
    config = (; nₘₐₓ, Ng)

    problem = LinearizationProblem(x₀; variables = vars) do x
        (; shape = Chebyshev(x[1], x[2], degree, complex(x[3], x[4])), λ = x[5])
    end

    function transition_vector_from_parameters(x)
        c = Chebyshev(x[1], x[2], degree, complex(x[3], x[4]))
        𝐓 = TransitionMatrices.transition_matrix(c, x[5], nₘₐₓ, Ng)
        chunks = [vcat(real.(vec(block)), imag.(vec(block))) for block in 𝐓.𝐓]
        return reduce(vcat, chunks)
    end

    function transition_vector_from_matrix(𝐓)
        chunks = [vcat(real.(vec(block)), imag.(vec(block))) for block in 𝐓.𝐓]
        return reduce(vcat, chunks)
    end

    support = supports_linearization(problem, EBCMLinearization(); config)
    @test Bool(support)

    result = linearize_transition_matrix(problem, EBCMLinearization(); config)
    @test result.metadata.backend == :ebcm_analytic
    @test !hasproperty(result.metadata, :reference)

    reference_jacobian = ForwardDiff.jacobian(transition_vector_from_parameters, x₀)
    for (j, var) in enumerate(vars)
        @test transition_vector_from_matrix(derivative(result, var))≈reference_jacobian[:, j] atol=1e-7 rtol=1e-7
    end
end

@testitem "EBCM cylinder linearization matches ForwardDiff references" begin
    using ForwardDiff

    nₘₐₓ = 3
    Ng = 28
    x₀ = [0.8, 1.4, 1.311, 0.02, 2π]
    vars = (:r, :h, :mᵣ, :mᵢ, :λ)
    config = (; nₘₐₓ, Ng)

    problem = LinearizationProblem(x₀; variables = vars) do x
        (; shape = Cylinder(x[1], x[2], complex(x[3], x[4])), λ = x[5])
    end

    function transition_vector_from_parameters(x)
        c = Cylinder(x[1], x[2], complex(x[3], x[4]))
        𝐓 = TransitionMatrices.transition_matrix(c, x[5], nₘₐₓ, Ng)
        chunks = [vcat(real.(vec(block)), imag.(vec(block))) for block in 𝐓.𝐓]
        return reduce(vcat, chunks)
    end

    function transition_vector_from_matrix(𝐓)
        chunks = [vcat(real.(vec(block)), imag.(vec(block))) for block in 𝐓.𝐓]
        return reduce(vcat, chunks)
    end

    support = supports_linearization(problem, EBCMLinearization(); config)
    @test Bool(support)

    result = linearize_transition_matrix(problem, EBCMLinearization(); config)
    @test result.metadata.backend == :ebcm_analytic
    @test !hasproperty(result.metadata, :reference)

    reference_jacobian = ForwardDiff.jacobian(transition_vector_from_parameters, x₀)
    for (j, var) in enumerate(vars)
        @test transition_vector_from_matrix(derivative(result, var))≈reference_jacobian[:, j] atol=1e-7 rtol=1e-7
    end
end

@testitem "IITM axisymmetric fixed-geometry linearization matches ForwardDiff references" begin
    using ForwardDiff

    nₘₐₓ = 2
    Nr = 5
    Nϑ = 12
    a = 0.9
    c = 1.2
    x₀ = [2π, 0.02, 1.311]
    vars = (:λ, :mᵢ, :mᵣ)
    config = (; nₘₐₓ, Nr, Nϑ)

    problem = LinearizationProblem(x₀; variables = vars) do x
        base = zero(x[1]) + zero(x[2]) + zero(x[3])
        (; shape = Spheroid(base + a, base + c, complex(x[3], x[2])),
            λ = x[1])
    end

    function transition_vector_from_matrix(𝐓)
        chunks = [vcat(real.(vec(block)), imag.(vec(block))) for block in 𝐓.𝐓]
        return reduce(vcat, chunks)
    end

    function transition_vector_from_parameters(x)
        base = zero(x[1]) + zero(x[2]) + zero(x[3])
        s = Spheroid(base + a, base + c, complex(x[3], x[2]))
        𝐓 = TransitionMatrices.transition_matrix_iitm(s, x[1], nₘₐₓ, Nr, Nϑ)
        return transition_vector_from_matrix(𝐓)
    end

    support = supports_linearization(problem, IITMLinearization(:axisymmetric); config)
    @test Bool(support)

    result = linearize_transition_matrix(problem, IITMLinearization(:axisymmetric);
        config)
    @test result.metadata.backend == :iitm_axisymmetric_analytic
    @test !hasproperty(result.metadata, :reference)

    reference = TransitionMatrices.transition_matrix_iitm(
        Spheroid(a, c,
            complex(x₀[3],
                x₀[2])),
        x₀[1], nₘₐₓ, Nr, Nϑ)
    reference_jacobian = ForwardDiff.jacobian(transition_vector_from_parameters, x₀)

    @test transition_vector_from_matrix(result.value) ≈
          transition_vector_from_matrix(reference)
    for (j, var) in enumerate(vars)
        @test transition_vector_from_matrix(derivative(result, var))≈reference_jacobian[:, j] atol=1e-7 rtol=1e-7
    end

    geometry_problem = LinearizationProblem([a]; variables = (:a,)) do x
        (; shape = Spheroid(x[1], c, complex(x₀[3], x₀[2])), λ = x₀[1])
    end
    geometry_support = supports_linearization(geometry_problem,
        IITMLinearization(:axisymmetric);
        config)
    @test !Bool(geometry_support)
    @test occursin("fixed-geometry", geometry_support.reason)
end

@testitem "IITM arbitrary fixed-geometry linearization matches ForwardDiff references" begin
    using ForwardDiff

    struct LinearizationArbitraryPrism{N, T, CT} <: AbstractShape{T, CT}
        s::Prism{N, T, CT}
        m::CT
    end

    TransitionMatrices.rmin(s::LinearizationArbitraryPrism) = rmin(s.s)
    TransitionMatrices.rmax(s::LinearizationArbitraryPrism) = rmax(s.s)
    TransitionMatrices.refractive_index(s::LinearizationArbitraryPrism, x) = refractive_index(s.s, x)

    nₘₐₓ = 2
    Nr = 3
    Nϑ = 8
    Nφ = 12
    x₀ = [2π, 1.311, 0.02]
    vars = (:λ, :mᵣ, :mᵢ)
    config = (; nₘₐₓ, Nr, Nϑ, Nφ)

    function shape_from_parameters(x)
        base = zero(x[1]) + zero(x[2]) + zero(x[3])
        prism = Prism(5, base + 0.8, base + 1.1, complex(x[2], x[3]))
        return LinearizationArbitraryPrism(prism, prism.m)
    end

    problem = LinearizationProblem(x₀; variables = vars) do x
        (; shape = shape_from_parameters(x), λ = x[1])
    end

    transition_vector_from_matrix(𝐓) = vcat(real.(vec(𝐓.container)), imag.(vec(𝐓.container)))

    function transition_vector_from_parameters(x)
        𝐓 = TransitionMatrices.transition_matrix_iitm(shape_from_parameters(x), x[1],
            nₘₐₓ, Nr, Nϑ, Nφ)
        return transition_vector_from_matrix(𝐓)
    end

    support = supports_linearization(problem, IITMLinearization(:arbitrary); config)
    @test Bool(support)

    result = linearize_transition_matrix(problem, IITMLinearization(:arbitrary);
        config)
    @test result.metadata.backend == :iitm_arbitrary_analytic
    @test !hasproperty(result.metadata, :reference)

    reference = TransitionMatrices.transition_matrix_iitm(shape_from_parameters(x₀),
        x₀[1], nₘₐₓ, Nr, Nϑ,
        Nφ)
    reference_jacobian = ForwardDiff.jacobian(transition_vector_from_parameters, x₀)

    @test transition_vector_from_matrix(result.value) ≈
          transition_vector_from_matrix(reference)
    for (j, var) in enumerate(vars)
        @test transition_vector_from_matrix(derivative(result, var))≈reference_jacobian[:, j] atol=1e-7 rtol=1e-7
    end
end

@testitem "IITM nfold fixed-geometry linearization matches production finite differences" begin
    nₘₐₓ = 2
    Nr = 3
    Nϑ = 8
    Nφ = 12
    x₀ = [2π, 0.02, 1.311]
    vars = (:λ, :mᵢ, :mᵣ)
    config = (; nₘₐₓ, Nr, Nϑ, Nφ)

    function shape_from_parameters(x)
        base = zero(x[1]) + zero(x[2]) + zero(x[3])
        return Prism(5, base + 0.8, base + 1.1, complex(x[3], x[2]))
    end

    problem = LinearizationProblem(x₀; variables = vars) do x
        (; shape = shape_from_parameters(x), λ = x[1])
    end

    transition_vector_from_matrix(𝐓) = vcat(real.(vec(𝐓.container)), imag.(vec(𝐓.container)))

    function transition_vector_from_parameters(x)
        𝐓 = TransitionMatrices.transition_matrix_iitm(shape_from_parameters(x), x[1],
            nₘₐₓ, Nr, Nϑ, Nφ)
        return transition_vector_from_matrix(𝐓)
    end

    support = supports_linearization(problem, IITMLinearization(:nfold); config)
    @test Bool(support)

    result = linearize_transition_matrix(problem, IITMLinearization(:nfold);
        config)
    @test result.metadata.backend == :iitm_nfold_analytic
    @test !hasproperty(result.metadata, :reference)

    reference = TransitionMatrices.transition_matrix_iitm(shape_from_parameters(x₀),
        x₀[1], nₘₐₓ, Nr, Nϑ,
        Nφ)
    @test transition_vector_from_matrix(result.value) ≈
          transition_vector_from_matrix(reference)

    δ = 1e-6
    for (j, var) in enumerate(vars)
        x⁺ = copy(x₀)
        x⁻ = copy(x₀)
        x⁺[j] += δ
        x⁻[j] -= δ
        ∂ref = (transition_vector_from_parameters(x⁺) .-
                transition_vector_from_parameters(x⁻)) ./ (2δ)
        @test transition_vector_from_matrix(derivative(result, var))≈∂ref atol=2e-5 rtol=2e-5
    end
end
