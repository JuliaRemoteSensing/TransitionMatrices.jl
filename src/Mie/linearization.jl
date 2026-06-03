const _MIE_LINEARIZATION_SUPPORTED_OUTPUTS = (
    :transition_matrix,
    :scattering_cross_section,
    :extinction_cross_section,
    :absorption_cross_section,
    :albedo,
    :amplitude_matrix,
)
const _MIE_LINEARIZATION_VARIABLES = (:x, :mᵣ, :mᵢ, :λ)

function _mie_linearization_input(problem::LinearizationProblem, config)
    rebuilt = rebuild(problem)
    x = _linearization_property(config, :x;
                                default = _linearization_property(rebuilt, :x))
    m = _linearization_property(config, :m;
                                default = _linearization_property(rebuilt, :m))
    nₘₐₓ = _linearization_property(config, :nₘₐₓ;
                                   default = _linearization_property(rebuilt, :nₘₐₓ))
    λ = _linearization_property(config, :λ;
                                default = _linearization_property(rebuilt, :λ;
                                                                  default = 2π))

    if isnothing(x) || isnothing(m) || isnothing(nₘₐₓ)
        return nothing
    end

    return (; x, m, nₘₐₓ = Int(nₘₐₓ), λ)
end

function _mie_linearization_variable_derivatives(problem::LinearizationProblem)
    vars = variables(problem)
    _linearization_variables_supported(vars, _MIE_LINEARIZATION_VARIABLES) || return nothing

    p = length(vars)
    ẋ = zeros(Float64, p)
    ṁ = zeros(ComplexF64, p)
    λ̇ = zeros(Float64, p)

    for (i, var) in enumerate(vars)
        if var == :x
            ẋ[i] = 1
        elseif var == :mᵣ
            ṁ[i] = 1
        elseif var == :mᵢ
            ṁ[i] = im
        elseif var == :λ
            λ̇[i] = 1
        end
    end

    return ẋ, ṁ, λ̇
end

function _mie_amplitude_angles(config)
    angles = _linearization_property(config, :angles)
    !isnothing(angles) && return Tuple(angles)

    ϑᵢ = _linearization_property(config, :ϑᵢ)
    φᵢ = _linearization_property(config, :φᵢ)
    ϑₛ = _linearization_property(config, :ϑₛ)
    φₛ = _linearization_property(config, :φₛ)
    any(isnothing, (ϑᵢ, φᵢ, ϑₛ, φₛ)) && return nothing
    return ϑᵢ, φᵢ, ϑₛ, φₛ
end

function supports_linearization(problem::LinearizationProblem, ::MieLinearization;
                                output::Symbol = :transition_matrix,
                                config = nothing)
    input = try
        _mie_linearization_input(problem, config)
    catch err
        return LinearizationSupport(false, "failed to rebuild Mie input: $err")
    end
    isnothing(input) &&
        return LinearizationSupport(false,
                                    "Mie linearization requires x, m, and nₘₐₓ")
    isnothing(_mie_linearization_variable_derivatives(problem)) &&
        return LinearizationSupport(false,
                                    "Mie linearization supports unique canonical variables drawn from :x, :mᵣ, :mᵢ, and :λ")
    output in _MIE_LINEARIZATION_SUPPORTED_OUTPUTS ||
        return LinearizationSupport(false, "Mie linearization does not support output :$output")
    output == :amplitude_matrix && isnothing(_mie_amplitude_angles(config)) &&
        return LinearizationSupport(false,
                                    "Mie amplitude matrix linearization requires angle config")

    return LinearizationSupport(true, "")
end

function _bhmie_linearized(::Type{T}, x, m, nₘₐₓ::Integer, ẋ, ṁ) where {T}
    C = complex(T)
    x = T(x)
    m = C(m)
    p = length(ẋ)
    ẋ = T.(ẋ)
    ṁ = C.(ṁ)

    y = m * x
    ẏ = @. ṁ * x + m * ẋ
    d = zeros(C, nₘₐₓ + 15)
    ḋ = zeros(C, nₘₐₓ + 15, p)

    for n in (nₘₐₓ + 14):-1:1
        c = T(n + 1)
        u = c / y
        v = d[n + 1] + u
        d[n] = u - inv(v)

        for j in 1:p
            u̇ = -c * ẏ[j] / y^2
            v̇ = ḋ[n + 1, j] + u̇
            ḋ[n, j] = u̇ + v̇ / v^2
        end
    end

    ψ₀ = cos(x)
    ψ₁ = sin(x)
    χ₀ = -sin(x)
    χ₁ = cos(x)
    ξ₁ = complex(ψ₁, -χ₁)

    ψ̇₀ = @. -sin(x) * ẋ
    ψ̇₁ = @. cos(x) * ẋ
    χ̇₀ = @. -cos(x) * ẋ
    χ̇₁ = @. -sin(x) * ẋ
    ξ̇₁ = C[complex(ψ̇₁[j], -χ̇₁[j]) for j in 1:p]
    ψ̇_next = similar(ψ̇₁)
    χ̇_next = similar(χ̇₁)

    a = zeros(C, nₘₐₓ)
    b = zeros(C, nₘₐₓ)
    ȧ = zeros(C, nₘₐₓ, p)
    ḃ = zeros(C, nₘₐₓ, p)

    for n in 1:nₘₐₓ
        c = T(2n - 1)
        ψ = c * ψ₁ / x - ψ₀
        χ = c * χ₁ / x - χ₀
        ξ = complex(ψ, -χ)
        A = d[n] / m + n / x
        B = d[n] * m + n / x
        numerator_a = A * ψ - ψ₁
        denominator_a = A * ξ - ξ₁
        numerator_b = B * ψ - ψ₁
        denominator_b = B * ξ - ξ₁

        for j in 1:p
            ψ̇ = c * (ψ̇₁[j] * x - ψ₁ * ẋ[j]) / x^2 - ψ̇₀[j]
            χ̇ = c * (χ̇₁[j] * x - χ₁ * ẋ[j]) / x^2 - χ̇₀[j]
            ξ̇ = complex(ψ̇, -χ̇)
            ψ̇_next[j] = ψ̇
            χ̇_next[j] = χ̇

            Ȧ = ḋ[n, j] / m - d[n] * ṁ[j] / m^2 - n * ẋ[j] / x^2
            numerator_dot = Ȧ * ψ + A * ψ̇ - ψ̇₁[j]
            denominator_dot = Ȧ * ξ + A * ξ̇ - ξ̇₁[j]
            ȧ[n, j] = (numerator_dot * denominator_a -
                       numerator_a * denominator_dot) / denominator_a^2

            Ḃ = ḋ[n, j] * m + d[n] * ṁ[j] - n * ẋ[j] / x^2
            numerator_dot = Ḃ * ψ + B * ψ̇ - ψ̇₁[j]
            denominator_dot = Ḃ * ξ + B * ξ̇ - ξ̇₁[j]
            ḃ[n, j] = (numerator_dot * denominator_b -
                       numerator_b * denominator_dot) / denominator_b^2
        end

        a[n] = numerator_a / denominator_a
        b[n] = numerator_b / denominator_b

        ψ₀, ψ₁ = ψ₁, ψ
        χ₀, χ₁ = χ₁, χ
        ξ₁ = complex(ψ₁, -χ₁)
        for j in 1:p
            ψ̇₀[j], ψ̇₁[j] = ψ̇₁[j], ψ̇_next[j]
            χ̇₀[j], χ̇₁[j] = χ̇₁[j], χ̇_next[j]
            ξ̇₁[j] = complex(ψ̇₁[j], -χ̇₁[j])
        end
    end

    return a, b, ȧ, ḃ
end

function _checked_mie_linearization_input(problem::LinearizationProblem,
                                          backend::MieLinearization,
                                          output::Symbol,
                                          config)
    output in _MIE_LINEARIZATION_SUPPORTED_OUTPUTS ||
        throw(UnsupportedLinearization(backend, output,
                                       "Mie linearization does not support output :$output"))

    input = try
        _mie_linearization_input(problem, config)
    catch err
        throw(UnsupportedLinearization(backend, output,
                                       "failed to rebuild Mie input: $err"))
    end
    isnothing(input) &&
        throw(UnsupportedLinearization(backend, output,
                                       "Mie linearization requires x, m, and nₘₐₓ"))

    variable_derivatives = _mie_linearization_variable_derivatives(problem)
    isnothing(variable_derivatives) &&
        throw(UnsupportedLinearization(backend, output,
                                       "Mie linearization supports unique canonical variables drawn from :x, :mᵣ, :mᵢ, and :λ"))
    output == :amplitude_matrix && isnothing(_mie_amplitude_angles(config)) &&
        throw(UnsupportedLinearization(backend, output,
                                       "Mie amplitude matrix linearization requires angle config"))

    return input, variable_derivatives
end

function _mie_transition_linearization_data(problem::LinearizationProblem,
                                            backend::MieLinearization,
                                            config;
                                            output::Symbol = :transition_matrix)
    input, (ẋ, ṁ, λ̇) = _checked_mie_linearization_input(problem, backend, output, config)
    T = typeof(float(real(input.x)))
    CT = complex(T)
    N = input.nₘₐₓ
    a, b, ȧ, ḃ = _bhmie_linearized(T, input.x, input.m, N, ẋ, ṁ)

    value = MieTransitionMatrix{CT, N, Vector{CT}}(a, b)
    jacobian = [MieTransitionMatrix{CT, N, Vector{CT}}(ȧ[:, j], ḃ[:, j])
                for j in 1:length(variables(problem))]

    result = LinearizationResult(value, jacobian, variables(problem);
                                 metadata = (; backend = :mie, λ = input.λ))
    return result, input, λ̇
end

function _mie_transition_linearization(problem::LinearizationProblem,
                                       backend::MieLinearization, config)
    result, _, _ = _mie_transition_linearization_data(problem, backend, config)
    return result
end

function linearize_transition_matrix(problem::LinearizationProblem,
                                     backend::MieLinearization; config = nothing)
    return _mie_transition_linearization(problem, backend, config)
end

function _mie_cross_section_linearization(result::LinearizationResult, λ, λ̇, kind::Symbol)
    mie = result.value
    vars = variables(result)
    N = order(mie)
    p = length(vars)
    λ = real(λ)
    scale = λ^2 / 2π
    scalė = @. λ * λ̇ / π
    base = zero(typeof(scale))
    basė = zeros(typeof(scale), p)

    for n in 1:N
        weight = 2n + 1
        if kind == :scattering
            base += weight * (abs2(mie.a[n]) + abs2(mie.b[n]))
        else
            base += weight * real(mie.a[n] + mie.b[n])
        end

        for j in 1:p
            ∂T = result.jacobian[j]
            if kind == :scattering
                basė[j] += weight * 2real(conj(mie.a[n]) * ∂T.a[n] +
                                           conj(mie.b[n]) * ∂T.b[n])
            else
                basė[j] += weight * real(∂T.a[n] + ∂T.b[n])
            end
        end
    end

    jacobian = @. basė * scale + base * scalė
    return base * scale, jacobian
end

function _mie_cross_section_pair_linearization(result::LinearizationResult, λ, λ̇)
    mie = result.value
    vars = variables(result)
    N = order(mie)
    p = length(vars)
    λ = real(λ)
    scale = λ^2 / 2π
    scalė = @. λ * λ̇ / π
    scattering_base = zero(typeof(scale))
    extinction_base = zero(typeof(scale))
    scattering_basė = zeros(typeof(scale), p)
    extinction_basė = zeros(typeof(scale), p)

    for n in 1:N
        weight = 2n + 1
        scattering_base += weight * (abs2(mie.a[n]) + abs2(mie.b[n]))
        extinction_base += weight * real(mie.a[n] + mie.b[n])

        for j in 1:p
            ∂T = result.jacobian[j]
            scattering_basė[j] += weight * 2real(conj(mie.a[n]) * ∂T.a[n] +
                                                  conj(mie.b[n]) * ∂T.b[n])
            extinction_basė[j] += weight * real(∂T.a[n] + ∂T.b[n])
        end
    end

    scattering_jacobian = @. scattering_basė * scale + scattering_base * scalė
    extinction_jacobian = @. extinction_basė * scale + extinction_base * scalė
    return (; scattering = scattering_base * scale,
            scattering_jacobian,
            extinction = extinction_base * scale,
            extinction_jacobian)
end

function linearize_observable(::typeof(scattering_cross_section),
                              problem::LinearizationProblem,
                              backend::MieLinearization; config = nothing)
    result, input, λ̇ = _mie_transition_linearization_data(problem, backend, config;
                                                          output = :scattering_cross_section)
    value, jacobian = _mie_cross_section_linearization(result, input.λ, λ̇, :scattering)
    return LinearizationResult(value, jacobian, variables(problem);
                               metadata = (; backend = :mie,
                                           observable = :scattering_cross_section))
end

function linearize_observable(::typeof(extinction_cross_section),
                              problem::LinearizationProblem,
                              backend::MieLinearization; config = nothing)
    result, input, λ̇ = _mie_transition_linearization_data(problem, backend, config;
                                                          output = :extinction_cross_section)
    value, jacobian = _mie_cross_section_linearization(result, input.λ, λ̇, :extinction)
    return LinearizationResult(value, jacobian, variables(problem);
                               metadata = (; backend = :mie,
                                           observable = :extinction_cross_section))
end

function linearize_observable(::typeof(absorption_cross_section),
                              problem::LinearizationProblem,
                              backend::MieLinearization; config = nothing)
    result, input, λ̇ = _mie_transition_linearization_data(problem, backend, config;
                                                          output = :absorption_cross_section)
    cross_sections = _mie_cross_section_pair_linearization(result, input.λ, λ̇)
    return LinearizationResult(cross_sections.extinction - cross_sections.scattering,
                               cross_sections.extinction_jacobian -
                               cross_sections.scattering_jacobian,
                               variables(problem);
                               metadata = (; backend = :mie,
                                           observable = :absorption_cross_section))
end

function linearize_observable(::typeof(albedo), problem::LinearizationProblem,
                              backend::MieLinearization; config = nothing)
    result, input, λ̇ = _mie_transition_linearization_data(problem, backend, config;
                                                          output = :albedo)
    cross_sections = _mie_cross_section_pair_linearization(result, input.λ, λ̇)
    value = cross_sections.scattering / cross_sections.extinction
    jacobian = @. (cross_sections.scattering_jacobian * cross_sections.extinction -
                  cross_sections.scattering * cross_sections.extinction_jacobian) /
                 cross_sections.extinction^2
    return LinearizationResult(value, jacobian, variables(problem);
                               metadata = (; backend = :mie, observable = :albedo))
end

function _mie_amplitude_matrix_linearization(result::LinearizationResult, λ, λ̇, angles)
    mie = result.value
    CT = eltype(mie.a)
    T = real(CT)
    N = order(mie)
    p = length(variables(result))
    scale = real(λ) / 2π
    𝐒₁₁, 𝐒₁₂, 𝐒₂₁, 𝐒₂₂ = zero(CT), zero(CT), zero(CT), zero(CT)
    raw_jacobian = zeros(CT, 2, 2, p)

    ϑᵢ, φᵢ, ϑₛ, φₛ = T.(angles)
    Δφ = φₛ - φᵢ

    πᵢ = OffsetArray(zeros(T, 2N + 1, N + 1), (-N):N, 0:N)
    τᵢ = OffsetArray(zeros(T, 2N + 1, N + 1), (-N):N, 0:N)
    πₛ = OffsetArray(zeros(T, 2N + 1, N + 1), (-N):N, 0:N)
    τₛ = OffsetArray(zeros(T, 2N + 1, N + 1), (-N):N, 0:N)

    for m in 0:N
        wigner_d_recursion!(view(πᵢ, m, m:N), 0, m, N, ϑᵢ;
                            deriv = view(τᵢ, m, m:N))
        wigner_d_recursion!(view(πₛ, m, m:N), 0, m, N, ϑₛ;
                            deriv = view(τₛ, m, m:N))
    end

    for m in 0:N
        for n in m:N
            if n > 0
                πᵢ[m, n] = pi_func(T, m, n, ϑᵢ; d = πᵢ[m, n])
                πₛ[m, n] = pi_func(T, m, n, ϑₛ; d = πₛ[m, n])
            end
            if m > 0
                π_sign = iseven(m + 1) ? one(T) : -one(T)
                τ_sign = iseven(m) ? one(T) : -one(T)
                πᵢ[-m, n] = π_sign * πᵢ[m, n]
                πₛ[-m, n] = π_sign * πₛ[m, n]
                τᵢ[-m, n] = τ_sign * τᵢ[m, n]
                τₛ[-m, n] = τ_sign * τₛ[m, n]
            end
        end
    end

    @inbounds for n in 1:N
        T₁₁ = -mie.b[n]
        T₂₂ = -mie.a[n]
        αₙ = CT(-1im) * T(2n + 1) / T(n * (n + 1))

        for m in (-n):n
            expiφ = cis(m * Δφ)
            πᵢₘₙ = πᵢ[m, n]
            τᵢₘₙ = τᵢ[m, n]
            πₛₘₙ = πₛ[m, n]
            τₛₘₙ = τₛ[m, n]
            ππ = πₛₘₙ * πᵢₘₙ
            πτ = πₛₘₙ * τᵢₘₙ
            τπ = τₛₘₙ * πᵢₘₙ
            ττ = τₛₘₙ * τᵢₘₙ
            αexpiφ = αₙ * expiφ

            𝐒₁₁ += αexpiφ * (T₁₁ * ππ + T₂₂ * ττ)
            𝐒₁₂ += αexpiφ * (T₁₁ * πτ + T₂₂ * τπ)
            𝐒₂₁ += αexpiφ * (T₁₁ * τπ + T₂₂ * πτ)
            𝐒₂₂ += αexpiφ * (T₁₁ * ττ + T₂₂ * ππ)

            for j in 1:p
                ∂T = result.jacobian[j]
                ∂T₁₁ = -∂T.b[n]
                ∂T₂₂ = -∂T.a[n]
                raw_jacobian[1, 1, j] += αexpiφ * (∂T₁₁ * ππ + ∂T₂₂ * ττ)
                raw_jacobian[1, 2, j] += αexpiφ * (∂T₁₁ * πτ + ∂T₂₂ * τπ)
                raw_jacobian[2, 1, j] += αexpiφ * (∂T₁₁ * τπ + ∂T₂₂ * πτ)
                raw_jacobian[2, 2, j] += αexpiφ * (∂T₁₁ * ττ + ∂T₂₂ * ππ)
            end
        end
    end

    value = (@SMatrix [𝐒₁₁ 𝐒₁₂ / 1im; 𝐒₂₁ * 1im 𝐒₂₂]) * scale
    jacobian = zeros(CT, 2, 2, p)
    for j in 1:p
        jacobian[1, 1, j] = raw_jacobian[1, 1, j] * scale
        jacobian[1, 2, j] = raw_jacobian[1, 2, j] / 1im * scale
        jacobian[2, 1, j] = raw_jacobian[2, 1, j] * 1im * scale
        jacobian[2, 2, j] = raw_jacobian[2, 2, j] * scale
        if !iszero(λ̇[j])
            jacobian[:, :, j] .+= value * (λ̇[j] / real(λ))
        end
    end

    return value, jacobian
end

function linearize_observable(::typeof(amplitude_matrix),
                              problem::LinearizationProblem,
                              backend::MieLinearization; config = nothing)
    result, input, λ̇ = _mie_transition_linearization_data(problem, backend, config;
                                                          output = :amplitude_matrix)
    angles = _mie_amplitude_angles(config)
    value, jacobian = _mie_amplitude_matrix_linearization(result, input.λ, λ̇, angles)

    return LinearizationResult(value, jacobian, variables(problem);
                               metadata = (; backend = :mie,
                                           observable = :amplitude_matrix))
end
