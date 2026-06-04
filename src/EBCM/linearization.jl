const _EBCM_SPHEROID_LINEARIZATION_VARIABLES = (:a, :c, :mᵣ, :mᵢ, :λ)
const _EBCM_CHEBYSHEV_LINEARIZATION_VARIABLES = (:r₀, :ε, :mᵣ, :mᵢ, :λ)
const _EBCM_CYLINDER_LINEARIZATION_VARIABLES = (:r, :h, :mᵣ, :mᵢ, :λ)

function _ebcm_linearization_shape(shape, λ)
    return shape
end

function _ebcm_linearization_shape(shape::Spheroid, λ)
    base = zero(shape.a) + zero(shape.c) + zero(real(λ)) +
           zero(real(shape.m)) + zero(imag(shape.m))
    a = base + shape.a
    c = base + shape.c
    m = complex(base + real(shape.m), base + imag(shape.m))
    return Spheroid(a, c, m)
end

function _ebcm_linearization_shape(shape::Chebyshev, λ)
    base = zero(shape.r₀) + zero(shape.ε) + zero(real(λ)) +
           zero(real(shape.m)) + zero(imag(shape.m))
    r₀ = base + shape.r₀
    ε = base + shape.ε
    m = complex(base + real(shape.m), base + imag(shape.m))
    return Chebyshev(r₀, ε, shape.n, m)
end

function _ebcm_linearization_shape(shape::Cylinder, λ)
    base = zero(shape.r) + zero(shape.h) + zero(real(λ)) +
           zero(real(shape.m)) + zero(imag(shape.m))
    r = base + shape.r
    h = base + shape.h
    m = complex(base + real(shape.m), base + imag(shape.m))
    return Cylinder(r, h, m)
end

function _ebcm_linearization_input(problem::LinearizationProblem, config,
                                   x = problem.x)
    rebuilt = rebuild(problem, x)
    shape = _linearization_property(config, :shape;
                                    default = _linearization_property(rebuilt, :shape))
    λ = _linearization_property(config, :λ;
                                default = _linearization_property(rebuilt, :λ))
    nₘₐₓ = _linearization_property(config, :nₘₐₓ;
                                   default = _linearization_property(rebuilt, :nₘₐₓ))
    Ng = _linearization_property(config, :Ng;
                                 default = _linearization_property(rebuilt, :Ng))

    if isnothing(shape) || isnothing(λ) || isnothing(nₘₐₓ) || isnothing(Ng)
        return nothing
    end

    shape = _ebcm_linearization_shape(shape, λ)
    return (; shape, λ, nₘₐₓ = Int(nₘₐₓ), Ng = Int(Ng))
end

function _ebcm_matrices(input)
    𝐏₀, 𝐔₀, cache = ebcm_matrices_m₀(input.shape, input.λ, input.nₘₐₓ, input.Ng;
                                      reuse = true)
    CT = eltype(𝐏₀)
    𝐏s = Vector{Matrix{CT}}(undef, input.nₘₐₓ + 1)
    𝐔s = Vector{Matrix{CT}}(undef, input.nₘₐₓ + 1)
    𝐏s[1] = 𝐏₀
    𝐔s[1] = 𝐔₀

    for m in 1:input.nₘₐₓ
        𝐏s[m + 1], 𝐔s[m + 1] = ebcm_matrices_m(m, input.shape, input.λ,
                                                 input.nₘₐₓ, input.Ng;
                                                 cache)
    end

    return 𝐏s, 𝐔s
end

struct _EBCMLinearizedScalar{V, D}
    value::V
    derivative::D
end

_ebcm_linearized(value, derivative) = _EBCMLinearizedScalar(value, derivative)
_ebcm_linearized_zero(::Type{T}) where {T} = _ebcm_linearized(zero(T), zero(T))
_ebcm_linearized_entry(value, derivative, i, n) =
    _ebcm_linearized(value[i, n], derivative[i, n])

Base.:+(x::_EBCMLinearizedScalar, y::_EBCMLinearizedScalar) =
    _ebcm_linearized(x.value + y.value, x.derivative + y.derivative)
Base.:+(x::_EBCMLinearizedScalar, y) =
    _ebcm_linearized(x.value + y, x.derivative)
Base.:+(x, y::_EBCMLinearizedScalar) =
    _ebcm_linearized(x + y.value, y.derivative)

Base.:-(x::_EBCMLinearizedScalar) =
    _ebcm_linearized(-x.value, -x.derivative)
Base.:-(x::_EBCMLinearizedScalar, y::_EBCMLinearizedScalar) =
    _ebcm_linearized(x.value - y.value, x.derivative - y.derivative)
Base.:-(x::_EBCMLinearizedScalar, y) =
    _ebcm_linearized(x.value - y, x.derivative)
Base.:-(x, y::_EBCMLinearizedScalar) =
    _ebcm_linearized(x - y.value, -y.derivative)

Base.:*(x::_EBCMLinearizedScalar, y::_EBCMLinearizedScalar) =
    _ebcm_linearized(x.value * y.value,
                     x.derivative * y.value + x.value * y.derivative)
Base.:*(x::_EBCMLinearizedScalar, y) =
    _ebcm_linearized(x.value * y, x.derivative * y)
Base.:*(x, y::_EBCMLinearizedScalar) =
    _ebcm_linearized(x * y.value, x * y.derivative)

Base.:/(x::_EBCMLinearizedScalar, y::_EBCMLinearizedScalar) =
    _ebcm_linearized(x.value / y.value,
                     (x.derivative * y.value - x.value * y.derivative) / y.value^2)
Base.:/(x::_EBCMLinearizedScalar, y) =
    _ebcm_linearized(x.value / y, x.derivative / y)
Base.:/(x, y::_EBCMLinearizedScalar) =
    _ebcm_linearized(x / y.value, -x * y.derivative / y.value^2)

function Base.:^(x::_EBCMLinearizedScalar, n::Integer)
    n == 0 && return _ebcm_linearized(one(x.value), zero(x.derivative))
    return _ebcm_linearized(x.value^n, n * x.value^(n - 1) * x.derivative)
end

function _ebcm_spheroid_radius_derivatives(s::Spheroid, x, r, r′, variable::Symbol)
    T = typeof(zero(s.a) + zero(s.c))
    ∂r = zeros(T, length(r))
    ∂r′ = zeros(T, length(r′))
    variable in (:a, :c) || return ∂r, ∂r′

    a = s.a
    c = s.c
    for i in eachindex(x)
        cosϑ = x[i]
        sinϑ = √(one(T) - cosϑ^2)
        denominator = a^2 * cosϑ^2 + c^2 * sinϑ^2
        common = cosϑ * sinϑ
        shape_factor = common * (inv(c^2) - inv(a^2))

        ∂r_a = r[i] * (inv(a) - a * cosϑ^2 / denominator)
        ∂r_c = r[i] * (inv(c) - c * sinϑ^2 / denominator)
        ∂shape_factor_a = 2common / a^3
        ∂shape_factor_c = -2common / c^3
        ∂r′_a = 3r[i]^2 * ∂r_a * shape_factor + r[i]^3 * ∂shape_factor_a
        ∂r′_c = 3r[i]^2 * ∂r_c * shape_factor + r[i]^3 * ∂shape_factor_c

        if variable == :a
            ∂r[i] = ∂r_a
            ∂r′[i] = ∂r′_a
        else
            ∂r[i] = ∂r_c
            ∂r′[i] = ∂r′_c
        end
    end

    return ∂r, ∂r′
end

function _ebcm_radius_derivatives(s::Spheroid, x, r, r′, variable::Symbol)
    return _ebcm_spheroid_radius_derivatives(s, x, r, r′, variable)
end

function _ebcm_radius_derivatives(c::Chebyshev, x, r, r′, variable::Symbol)
    T = typeof(zero(c.r₀) + zero(c.ε))
    ∂r = zeros(T, length(r))
    ∂r′ = zeros(T, length(r′))
    variable in (:r₀, :ε) || return ∂r, ∂r′

    for i in eachindex(x)
        nϑ = acos(x[i]) * c.n
        cosnϑ = cos(nϑ)
        sinnϑ = sin(nϑ)

        if variable == :r₀
            ∂r[i] = one(T) + c.ε * cosnϑ
            ∂r′[i] = -c.ε * c.n * sinnϑ
        else
            ∂r[i] = c.r₀ * cosnϑ
            ∂r′[i] = -c.r₀ * c.n * sinnϑ
        end
    end

    return ∂r, ∂r′
end

function _ebcm_radius_derivatives(c::Cylinder, x, r, r′, variable::Symbol)
    T = typeof(zero(c.r) + zero(c.h))
    return zeros(T, length(r)), zeros(T, length(r′))
end

function _ebcm_geometric_derivatives(s, x, w, r, r′, variable::Symbol)
    T = eltype(r)
    ∂x = zeros(T, length(x))
    ∂w = zeros(T, length(w))
    ∂r, ∂r′ = _ebcm_radius_derivatives(s, x, r, r′, variable)
    ∂ϑ = zeros(T, length(x))
    return ∂x, ∂w, ∂r, ∂r′, ∂ϑ
end

function _ebcm_geometric_derivatives(c::Cylinder, x, w, r, r′, variable::Symbol)
    T = typeof(zero(c.r) + zero(c.h))
    ∂x = zeros(T, length(x))
    ∂w = zeros(T, length(w))
    ∂r = zeros(T, length(r))
    ∂r′ = zeros(T, length(r′))
    ∂ϑ = zeros(T, length(x))
    variable in (:r, :h) || return ∂x, ∂w, ∂r, ∂r′, ∂ϑ

    R = c.r
    H = c.h / 2
    ∂R = variable == :r ? one(T) : zero(T)
    ∂H = variable == :h ? one(T) / 2 : zero(T)
    L = hypot(R, H)
    ∂xx = -∂H / L + H * (R * ∂R + H * ∂H) / L^3

    ng = length(x) ÷ 2
    ng1 = ng ÷ 2
    ng2 = ng - ng1
    x1, w1 = gausslegendre(T, ng1)
    x2, w2 = gausslegendre(T, ng2)

    @. ∂x[1:ng1] = 0.5 * (x1 + 1) * ∂xx
    @. ∂w[1:ng1] = 0.5 * ∂xx * w1
    @. ∂x[(ng1 + 1):ng] = 0.5 * (1 - x2) * ∂xx
    @. ∂w[(ng1 + 1):ng] = -0.5 * ∂xx * w2

    for i in 1:ng
        j = length(x) + 1 - i
        ∂x[j] = -∂x[i]
        ∂w[j] = ∂w[i]

        cosϑ = abs(x[i])
        ∂cosϑ = x[i] < 0 ? -∂x[i] : ∂x[i]
        sinϑ = √(one(T) - cosϑ^2)
        ∂sinϑ = -cosϑ * ∂cosϑ / sinϑ

        if H / cosϑ < R / sinϑ
            ∂rᵢ = ∂H / cosϑ - H * ∂cosϑ / cosϑ^2
            ∂r′ᵣₐw = ∂H * sinϑ / cosϑ^2 +
                      H * (∂sinϑ / cosϑ^2 - 2sinϑ * ∂cosϑ / cosϑ^3)
        else
            ∂rᵢ = ∂R / sinϑ - R * ∂sinϑ / sinϑ^2
            ∂r′ᵣₐw = -∂R * cosϑ / sinϑ^2 - R * ∂cosϑ / sinϑ^2 +
                      2R * cosϑ * ∂sinϑ / sinϑ^3
        end

        ∂r[i] = ∂rᵢ
        ∂r[j] = ∂rᵢ
        ∂r′[i] = -∂r′ᵣₐw
        ∂r′[j] = ∂r′ᵣₐw
    end

    for i in eachindex(x)
        ∂ϑ[i] = -∂x[i] / √(one(T) - x[i]^2)
    end

    return ∂x, ∂w, ∂r, ∂r′, ∂ϑ
end

function _ebcm_shape_real_zero(s::Spheroid)
    return zero(s.a) + zero(s.c)
end

function _ebcm_shape_real_zero(c::Chebyshev)
    return zero(c.r₀) + zero(c.ε)
end

function _ebcm_shape_real_zero(c::Cylinder)
    return zero(c.r) + zero(c.h)
end

function _ebcm_parameter_derivatives(input, variable::Symbol)
    s = input.shape
    T = typeof(_ebcm_shape_real_zero(s) + zero(real(input.λ)))
    CT = typeof(complex(zero(T) + real(s.m), zero(T) + imag(s.m)))
    k = 2 * T(π) / input.λ
    ∂λ = variable == :λ ? one(T) : zero(T)
    ∂k = -k * ∂λ / input.λ
    ∂m = if variable == :mᵣ
        one(CT)
    elseif variable == :mᵢ
        one(T) * im
    else
        zero(CT)
    end
    return k, ∂k, ∂m
end

function _ebcm_ricatti_argument_derivatives!(∂f, ∂f′, f, f′, z, ∂z)
    for n in eachindex(f, f′)
        ∂f[n] = f′[n] * ∂z
        ∂f′[n] = (n * (n + 1) / z^2 - 1) * f[n] * ∂z
    end
end

function _ebcm_wigner_tables(::Type{T}, m::Integer, nₘₐₓ::Integer, Ng::Integer,
                             ϑ, nₘᵢₙ::Integer, ∂ϑ = nothing) where {T}
    d = OffsetArray(zeros(T, Ng, nₘₐₓ - m + 1), 1:Ng, m:nₘₐₓ)
    𝜋 = similar(d)
    τ = similar(d)
    ∂d = similar(d)
    ∂𝜋 = similar(d)
    ∂τ = similar(d)
    fill!(∂d, zero(T))
    fill!(∂𝜋, zero(T))
    fill!(∂τ, zero(T))

    for i in eachindex(ϑ)
        wigner_d_recursion!(view(d, i, :), 0, m, nₘₐₓ, ϑ[i];
                            deriv = view(τ, i, :))

        sinϑ = sin(ϑ[i])
        cosϑ = cos(ϑ[i])
        ∂ϑᵢ = isnothing(∂ϑ) ? zero(T) : ∂ϑ[i]
        for n in nₘᵢₙ:nₘₐₓ
            𝜋[i, n] = pi_func(T, m, n, ϑ[i]; d = d[i, n])
            ∂d[i, n] = τ[i, n] * ∂ϑᵢ
            d² = -cosϑ / sinϑ * τ[i, n] -
                 (n * (n + 1) - m^2 / sinϑ^2) * d[i, n]
            ∂τ[i, n] = d² * ∂ϑᵢ
            ∂𝜋[i, n] = m * (∂d[i, n] / sinϑ -
                           d[i, n] * cosϑ * ∂ϑᵢ / sinϑ^2)
        end
    end

    return d, 𝜋, τ, ∂d, ∂𝜋, ∂τ
end

function _ebcm_directional_data(input, cache, variable::Symbol)
    x, w, r, r′, ϑ, a, A, ψ, ψ′, χ, χ′, ψₛ, ψₛ′, _, _ = cache
    s = input.shape
    nₘₐₓ = input.nₘₐₓ
    ng = size(ψ, 1)
    CT = eltype(ψₛ)
    k, ∂k, ∂m = _ebcm_parameter_derivatives(input, variable)
    ∂x, ∂w, ∂r, ∂r′, ∂ϑ = _ebcm_geometric_derivatives(s, x, w, r, r′, variable)

    ∂ψ = zeros(CT, ng, nₘₐₓ)
    ∂ψ′ = similar(∂ψ)
    ∂χ = similar(∂ψ)
    ∂χ′ = similar(∂ψ)
    ∂ψₛ = similar(∂ψ)
    ∂ψₛ′ = similar(∂ψ)

    for i in 1:ng
        kr = k * r[i]
        ∂kr = ∂k * r[i] + k * ∂r[i]
        ksr = s.m * k * r[i]
        ∂ksr = ∂m * k * r[i] + s.m * ∂k * r[i] + s.m * k * ∂r[i]

        _ebcm_ricatti_argument_derivatives!(view(∂ψ, i, :),
                                            view(∂ψ′, i, :),
                                            view(ψ, i, :),
                                            view(ψ′, i, :), kr, ∂kr)
        _ebcm_ricatti_argument_derivatives!(view(∂χ, i, :),
                                            view(∂χ′, i, :),
                                            view(χ, i, :),
                                            view(χ′, i, :), kr, ∂kr)
        _ebcm_ricatti_argument_derivatives!(view(∂ψₛ, i, :),
                                            view(∂ψₛ′, i, :),
                                            view(ψₛ, i, :),
                                            view(ψₛ′, i, :), ksr, ∂ksr)
    end

    return (; x, w, r, r′, ϑ, a, A, ψ, ψ′, χ, χ′, ψₛ, ψₛ′,
            k, ∂k, ∂m, ∂r, ∂r′, ∂ψ, ∂ψ′,
            ∂χ, ∂χ′, ∂ψₛ, ∂ψₛ′, ∂x, ∂w, ∂ϑ)
end

function _ebcm_matrices_m₀_derivative(input, data)
    s = input.shape
    nₘₐₓ = input.nₘₐₓ
    Ng = input.Ng
    ng = size(data.ψ, 1)
    T = eltype(data.r)
    CT = eltype(data.ψₛ)
    d, 𝜋, τ, ∂d, ∂𝜋, ∂τ = _ebcm_wigner_tables(T, 0, nₘₐₓ, Ng, data.ϑ, 0,
                                                 data.∂ϑ)

    ∂𝐏 = zeros(CT, 2nₘₐₓ, 2nₘₐₓ)
    ∂𝐏₁₁ = view(∂𝐏, 1:nₘₐₓ, 1:nₘₐₓ)
    ∂𝐏₂₂ = view(∂𝐏, (nₘₐₓ + 1):(2nₘₐₓ), (nₘₐₓ + 1):(2nₘₐₓ))
    ∂𝐔 = zeros(CT, 2nₘₐₓ, 2nₘₐₓ)
    ∂𝐔₁₁ = view(∂𝐔, 1:nₘₐₓ, 1:nₘₐₓ)
    ∂𝐔₂₂ = view(∂𝐔, (nₘₐₓ + 1):(2nₘₐₓ), (nₘₐₓ + 1):(2nₘₐₓ))

    m = _ebcm_linearized(s.m, data.∂m)
    k = _ebcm_linearized(data.k, data.∂k)

    for n in 1:nₘₐₓ, n′ in 1:nₘₐₓ
        if isodd(n + n′)
            continue
        end

        if n != n′
            PL₁ = _ebcm_linearized_zero(CT)
            PL₂ = _ebcm_linearized_zero(CT)
            PL₇ = _ebcm_linearized_zero(CT)
            PL₈ = _ebcm_linearized_zero(CT)
            UL₁ = _ebcm_linearized_zero(CT)
            UL₂ = _ebcm_linearized_zero(CT)
            UL₇ = _ebcm_linearized_zero(CT)
            UL₈ = _ebcm_linearized_zero(CT)

            for i in 1:ng
                w = _ebcm_linearized(data.w[i], data.∂w[i])
                r = _ebcm_linearized(data.r[i], data.∂r[i])
                r′ = _ebcm_linearized(data.r′[i], data.∂r′[i])
                dn = _ebcm_linearized_entry(d, ∂d, i, n)
                dn′ = _ebcm_linearized_entry(d, ∂d, i, n′)
                τn = _ebcm_linearized_entry(τ, ∂τ, i, n)
                τn′ = _ebcm_linearized_entry(τ, ∂τ, i, n′)
                ψn = _ebcm_linearized(data.ψ[i, n], data.∂ψ[i, n])
                ψ′n = _ebcm_linearized(data.ψ′[i, n], data.∂ψ′[i, n])
                χn = _ebcm_linearized(data.χ[i, n], data.∂χ[i, n])
                χ′n = _ebcm_linearized(data.χ′[i, n], data.∂χ′[i, n])
                ψₛn′ = _ebcm_linearized(data.ψₛ[i, n′], data.∂ψₛ[i, n′])
                ψₛ′n′ = _ebcm_linearized(data.ψₛ′[i, n′], data.∂ψₛ′[i, n′])
                kr = k * r
                factor = w * k * r′

                PL₁ += factor * τn * dn′ * ψn * ψₛn′
                PL₂ += factor * dn * τn′ * ψn * ψₛn′
                PL₇ += factor * τn * dn′ *
                       (ψ′n * ψₛ′n′ +
                        data.a[n] * ψn * ψₛn′ / (m * kr^2))
                PL₈ += factor * dn * τn′ *
                       (ψ′n * ψₛ′n′ +
                        data.a[n′] * ψn * ψₛn′ / (m * kr^2))

                UL₁ += factor * τn * dn′ * χn * ψₛn′
                UL₂ += factor * dn * τn′ * χn * ψₛn′
                UL₇ += factor * τn * dn′ *
                       (χ′n * ψₛ′n′ +
                        data.a[n] * χn * ψₛn′ / (m * kr^2))
                UL₈ += factor * dn * τn′ *
                       (χ′n * ψₛ′n′ +
                        data.a[n′] * χn * ψₛn′ / (m * kr^2))
            end

            prefactor = 1im * data.A[n] * data.A[n′] * (m^2 - 1) /
                        (m * (data.a[n] - data.a[n′]))
            ∂𝐏₁₁[n, n′] = (prefactor * (data.a[n] * PL₂ -
                                        data.a[n′] * PL₁)).derivative
            ∂𝐏₂₂[n, n′] = (prefactor * (data.a[n] * PL₈ -
                                        data.a[n′] * PL₇)).derivative
            ∂𝐔₁₁[n, n′] = (prefactor * (data.a[n] * UL₂ -
                                        data.a[n′] * UL₁)).derivative
            ∂𝐔₂₂[n, n′] = (prefactor * (data.a[n] * UL₈ -
                                        data.a[n′] * UL₇)).derivative
        else
            PL̃₁ = _ebcm_linearized_zero(CT)
            PL̃₂ = _ebcm_linearized_zero(CT)
            PL̃₃ = _ebcm_linearized_zero(CT)
            UL̃₁ = _ebcm_linearized_zero(CT)
            UL̃₂ = _ebcm_linearized_zero(CT)
            UL̃₃ = _ebcm_linearized_zero(CT)

            for i in 1:ng
                w = _ebcm_linearized(data.w[i], data.∂w[i])
                r = _ebcm_linearized(data.r[i], data.∂r[i])
                r′ = _ebcm_linearized(data.r′[i], data.∂r′[i])
                dn = _ebcm_linearized_entry(d, ∂d, i, n)
                𝜋n = _ebcm_linearized_entry(𝜋, ∂𝜋, i, n)
                τn = _ebcm_linearized_entry(τ, ∂τ, i, n)
                ψn = _ebcm_linearized(data.ψ[i, n], data.∂ψ[i, n])
                ψ′n = _ebcm_linearized(data.ψ′[i, n], data.∂ψ′[i, n])
                χn = _ebcm_linearized(data.χ[i, n], data.∂χ[i, n])
                χ′n = _ebcm_linearized(data.χ′[i, n], data.∂χ′[i, n])
                ψₛn = _ebcm_linearized(data.ψₛ[i, n], data.∂ψₛ[i, n])
                ψₛ′n = _ebcm_linearized(data.ψₛ′[i, n], data.∂ψₛ′[i, n])
                kr = k * r

                PL̃₁ += w * (𝜋n^2 + τn^2) *
                        (ψ′n * ψₛn - m * ψn * ψₛ′n)
                PL̃₂ += w * (𝜋n^2 + τn^2) *
                        (m * ψ′n * ψₛn - ψn * ψₛ′n)
                PL̃₃ += w * k * r′ * τn * dn *
                        ψn * ψₛn / (m * kr^2)

                UL̃₁ += w * (𝜋n^2 + τn^2) *
                        (χ′n * ψₛn - m * χn * ψₛ′n)
                UL̃₂ += w * (𝜋n^2 + τn^2) *
                        (m * χ′n * ψₛn - χn * ψₛ′n)
                UL̃₃ += w * k * r′ * τn * dn *
                        χn * ψₛn / (m * kr^2)
            end

            ∂𝐏₁₁[n, n] = (-1im / m * data.A[n]^2 * PL̃₁).derivative
            ∂𝐏₂₂[n, n] = (-1im / m * data.A[n]^2 *
                           (PL̃₂ + (m^2 - 1) * data.a[n] * PL̃₃)).derivative
            ∂𝐔₁₁[n, n] = (-1im / m * data.A[n]^2 * UL̃₁).derivative
            ∂𝐔₂₂[n, n] = (-1im / m * data.A[n]^2 *
                           (UL̃₂ + (m^2 - 1) * data.a[n] * UL̃₃)).derivative
        end
    end

    return ∂𝐏, ∂𝐔
end

function _ebcm_matrices_m_derivative(input, data, m_order::Integer)
    s = input.shape
    nₘₐₓ = input.nₘₐₓ
    Ng = input.Ng
    nₘᵢₙ = max(1, m_order)
    nn = nₘₐₓ - nₘᵢₙ + 1
    ng = size(data.ψ, 1)
    T = eltype(data.r)
    CT = eltype(data.ψₛ)
    d, 𝜋, τ, ∂d, ∂𝜋, ∂τ = _ebcm_wigner_tables(T, m_order, nₘₐₓ, Ng,
                                                 data.ϑ, nₘᵢₙ, data.∂ϑ)

    ∂𝐏 = zeros(CT, 2nn, 2nn)
    ∂𝐏₁₁ = OffsetArray(view(∂𝐏, 1:nn, 1:nn), nₘᵢₙ:nₘₐₓ, nₘᵢₙ:nₘₐₓ)
    ∂𝐏₁₂ = OffsetArray(view(∂𝐏, 1:nn, (nn + 1):(2nn)), nₘᵢₙ:nₘₐₓ,
                        nₘᵢₙ:nₘₐₓ)
    ∂𝐏₂₁ = OffsetArray(view(∂𝐏, (nn + 1):(2nn), 1:nn), nₘᵢₙ:nₘₐₓ,
                        nₘᵢₙ:nₘₐₓ)
    ∂𝐏₂₂ = OffsetArray(view(∂𝐏, (nn + 1):(2nn), (nn + 1):(2nn)),
                        nₘᵢₙ:nₘₐₓ, nₘᵢₙ:nₘₐₓ)
    ∂𝐔 = zeros(CT, 2nn, 2nn)
    ∂𝐔₁₁ = OffsetArray(view(∂𝐔, 1:nn, 1:nn), nₘᵢₙ:nₘₐₓ, nₘᵢₙ:nₘₐₓ)
    ∂𝐔₁₂ = OffsetArray(view(∂𝐔, 1:nn, (nn + 1):(2nn)), nₘᵢₙ:nₘₐₓ,
                        nₘᵢₙ:nₘₐₓ)
    ∂𝐔₂₁ = OffsetArray(view(∂𝐔, (nn + 1):(2nn), 1:nn), nₘᵢₙ:nₘₐₓ,
                        nₘᵢₙ:nₘₐₓ)
    ∂𝐔₂₂ = OffsetArray(view(∂𝐔, (nn + 1):(2nn), (nn + 1):(2nn)),
                        nₘᵢₙ:nₘₐₓ, nₘᵢₙ:nₘₐₓ)

    m = _ebcm_linearized(s.m, data.∂m)
    k = _ebcm_linearized(data.k, data.∂k)

    for n in nₘᵢₙ:nₘₐₓ, n′ in nₘᵢₙ:nₘₐₓ
        if !iseven(n + n′)
            PK₁ = _ebcm_linearized_zero(CT)
            PK₂ = _ebcm_linearized_zero(CT)
            UK₁ = _ebcm_linearized_zero(CT)
            UK₂ = _ebcm_linearized_zero(CT)

            for i in 1:ng
                w = _ebcm_linearized(data.w[i], data.∂w[i])
                r′ = _ebcm_linearized(data.r′[i], data.∂r′[i])
                dn′ = _ebcm_linearized_entry(d, ∂d, i, n′)
                𝜋n = _ebcm_linearized_entry(𝜋, ∂𝜋, i, n)
                ψn = _ebcm_linearized(data.ψ[i, n], data.∂ψ[i, n])
                ψ′n = _ebcm_linearized(data.ψ′[i, n], data.∂ψ′[i, n])
                χn = _ebcm_linearized(data.χ[i, n], data.∂χ[i, n])
                χ′n = _ebcm_linearized(data.χ′[i, n], data.∂χ′[i, n])
                ψₛn′ = _ebcm_linearized(data.ψₛ[i, n′], data.∂ψₛ[i, n′])
                ψₛ′n′ = _ebcm_linearized(data.ψₛ′[i, n′], data.∂ψₛ′[i, n′])
                factor = w * k * r′ * 𝜋n * dn′

                PK₁ += factor * ψn * ψₛ′n′
                PK₂ += factor * ψ′n * ψₛn′
                UK₁ += factor * χn * ψₛ′n′
                UK₂ += factor * χ′n * ψₛn′
            end

            ∂𝐏₁₂[n, n′] = (data.A[n] * data.A[n′] * (m^2 - 1) / m * PK₁).derivative
            ∂𝐏₂₁[n, n′] = (data.A[n] * data.A[n′] * (1 - m^2) / m * PK₂).derivative
            ∂𝐔₁₂[n, n′] = (data.A[n] * data.A[n′] * (m^2 - 1) / m * UK₁).derivative
            ∂𝐔₂₁[n, n′] = (data.A[n] * data.A[n′] * (1 - m^2) / m * UK₂).derivative
        end

        if !isodd(n + n′)
            if n != n′
                PL₁ = _ebcm_linearized_zero(CT)
                PL₂ = _ebcm_linearized_zero(CT)
                PL₇ = _ebcm_linearized_zero(CT)
                PL₈ = _ebcm_linearized_zero(CT)
                UL₁ = _ebcm_linearized_zero(CT)
                UL₂ = _ebcm_linearized_zero(CT)
                UL₇ = _ebcm_linearized_zero(CT)
                UL₈ = _ebcm_linearized_zero(CT)

                for i in 1:ng
                    w = _ebcm_linearized(data.w[i], data.∂w[i])
                    r = _ebcm_linearized(data.r[i], data.∂r[i])
                    r′ = _ebcm_linearized(data.r′[i], data.∂r′[i])
                    dn = _ebcm_linearized_entry(d, ∂d, i, n)
                    dn′ = _ebcm_linearized_entry(d, ∂d, i, n′)
                    τn = _ebcm_linearized_entry(τ, ∂τ, i, n)
                    τn′ = _ebcm_linearized_entry(τ, ∂τ, i, n′)
                    ψn = _ebcm_linearized(data.ψ[i, n], data.∂ψ[i, n])
                    ψ′n = _ebcm_linearized(data.ψ′[i, n], data.∂ψ′[i, n])
                    χn = _ebcm_linearized(data.χ[i, n], data.∂χ[i, n])
                    χ′n = _ebcm_linearized(data.χ′[i, n], data.∂χ′[i, n])
                    ψₛn′ = _ebcm_linearized(data.ψₛ[i, n′], data.∂ψₛ[i, n′])
                    ψₛ′n′ = _ebcm_linearized(data.ψₛ′[i, n′], data.∂ψₛ′[i, n′])
                    kr = k * r
                    factor = w * k * r′

                    PL₁ += factor * τn * dn′ * ψn * ψₛn′
                    PL₂ += factor * dn * τn′ * ψn * ψₛn′
                    PL₇ += factor * τn * dn′ *
                           (ψ′n * ψₛ′n′ +
                            n * (n + 1) * ψn * ψₛn′ / (m * kr^2))
                    PL₈ += factor * dn * τn′ *
                           (ψ′n * ψₛ′n′ +
                            n′ * (n′ + 1) * ψn * ψₛn′ / (m * kr^2))

                    UL₁ += factor * τn * dn′ * χn * ψₛn′
                    UL₂ += factor * dn * τn′ * χn * ψₛn′
                    UL₇ += factor * τn * dn′ *
                           (χ′n * ψₛ′n′ +
                            n * (n + 1) * χn * ψₛn′ / (m * kr^2))
                    UL₈ += factor * dn * τn′ *
                           (χ′n * ψₛ′n′ +
                            n′ * (n′ + 1) * χn * ψₛn′ / (m * kr^2))
                end

                prefactor = 1im * data.A[n] * data.A[n′] * (m^2 - 1) /
                            (m * (data.a[n] - data.a[n′]))
                ∂𝐏₁₁[n, n′] = (prefactor * (data.a[n] * PL₂ -
                                            data.a[n′] * PL₁)).derivative
                ∂𝐏₂₂[n, n′] = (prefactor * (data.a[n] * PL₈ -
                                            data.a[n′] * PL₇)).derivative
                ∂𝐔₁₁[n, n′] = (prefactor * (data.a[n] * UL₂ -
                                            data.a[n′] * UL₁)).derivative
                ∂𝐔₂₂[n, n′] = (prefactor * (data.a[n] * UL₈ -
                                            data.a[n′] * UL₇)).derivative
            else
                PL̃₁ = _ebcm_linearized_zero(CT)
                PL̃₂ = _ebcm_linearized_zero(CT)
                PL̃₃ = _ebcm_linearized_zero(CT)
                UL̃₁ = _ebcm_linearized_zero(CT)
                UL̃₂ = _ebcm_linearized_zero(CT)
                UL̃₃ = _ebcm_linearized_zero(CT)

                for i in 1:ng
                    w = _ebcm_linearized(data.w[i], data.∂w[i])
                    r = _ebcm_linearized(data.r[i], data.∂r[i])
                    r′ = _ebcm_linearized(data.r′[i], data.∂r′[i])
                    dn = _ebcm_linearized_entry(d, ∂d, i, n)
                    𝜋n = _ebcm_linearized_entry(𝜋, ∂𝜋, i, n)
                    τn = _ebcm_linearized_entry(τ, ∂τ, i, n)
                    ψn = _ebcm_linearized(data.ψ[i, n], data.∂ψ[i, n])
                    ψ′n = _ebcm_linearized(data.ψ′[i, n], data.∂ψ′[i, n])
                    χn = _ebcm_linearized(data.χ[i, n], data.∂χ[i, n])
                    χ′n = _ebcm_linearized(data.χ′[i, n], data.∂χ′[i, n])
                    ψₛn = _ebcm_linearized(data.ψₛ[i, n], data.∂ψₛ[i, n])
                    ψₛ′n = _ebcm_linearized(data.ψₛ′[i, n], data.∂ψₛ′[i, n])
                    kr = k * r

                    PL̃₁ += w * (𝜋n^2 + τn^2) *
                            (ψ′n * ψₛn - m * ψn * ψₛ′n)
                    PL̃₂ += w * (𝜋n^2 + τn^2) *
                            (m * ψ′n * ψₛn - ψn * ψₛ′n)
                    PL̃₃ += w * k * r′ * τn * dn *
                            ψn * ψₛn / (m * kr^2)

                    UL̃₁ += w * (𝜋n^2 + τn^2) *
                            (χ′n * ψₛn - m * χn * ψₛ′n)
                    UL̃₂ += w * (𝜋n^2 + τn^2) *
                            (m * χ′n * ψₛn - χn * ψₛ′n)
                    UL̃₃ += w * k * r′ * τn * dn *
                            χn * ψₛn / (m * kr^2)
                end

                ∂𝐏₁₁[n, n] = (-1im / m * data.A[n]^2 * PL̃₁).derivative
                ∂𝐏₂₂[n, n] = (-1im / m * data.A[n]^2 *
                               (PL̃₂ + (m^2 - 1) * data.a[n] * PL̃₃)).derivative
                ∂𝐔₁₁[n, n] = (-1im / m * data.A[n]^2 * UL̃₁).derivative
                ∂𝐔₂₂[n, n] = (-1im / m * data.A[n]^2 *
                               (UL̃₂ + (m^2 - 1) * data.a[n] * UL̃₃)).derivative
            end
        end
    end

    return ∂𝐏, ∂𝐔
end

function _ebcm_matrix_derivatives(input, variables)
    _, _, cache = ebcm_matrices_m₀(input.shape, input.λ, input.nₘₐₓ, input.Ng;
                                   reuse = true)
    map(variables) do variable
        data = _ebcm_directional_data(input, cache, variable)
        ∂𝐏s = Vector{Matrix{eltype(data.ψₛ)}}(undef, input.nₘₐₓ + 1)
        ∂𝐔s = Vector{Matrix{eltype(data.ψₛ)}}(undef, input.nₘₐₓ + 1)
        ∂𝐏s[1], ∂𝐔s[1] = _ebcm_matrices_m₀_derivative(input, data)
        for m in 1:input.nₘₐₓ
            ∂𝐏s[m + 1], ∂𝐔s[m + 1] = _ebcm_matrices_m_derivative(input, data, m)
        end
        ∂𝐏s, ∂𝐔s
    end
end

function _ebcm_linearization_variables(::Spheroid)
    return _EBCM_SPHEROID_LINEARIZATION_VARIABLES
end

function _ebcm_linearization_variables(::Chebyshev)
    return _EBCM_CHEBYSHEV_LINEARIZATION_VARIABLES
end

function _ebcm_linearization_variables(::Cylinder)
    return _EBCM_CYLINDER_LINEARIZATION_VARIABLES
end

function _ebcm_linearization_variables(shape)
    return nothing
end

function _ebcm_variable_list_message(canonical)
    return join((":" * String(var) for var in canonical), ", ", " and ")
end

function supports_linearization(problem::LinearizationProblem, ::EBCMLinearization;
                                output::Symbol = :transition_matrix,
                                config = nothing)
    output == :transition_matrix ||
        return LinearizationSupport(false,
                                    "EBCM analytical linearization only supports transition matrices")

    input = try
        _ebcm_linearization_input(problem, config)
    catch err
        return LinearizationSupport(false, "failed to rebuild EBCM input: $err")
    end
    isnothing(input) &&
        return LinearizationSupport(false,
                                    "EBCM analytical linearization requires shape, λ, nₘₐₓ, and Ng")
    canonical_variables = _ebcm_linearization_variables(input.shape)
    isnothing(canonical_variables) &&
        return LinearizationSupport(false,
                                    "EBCM analytical linearization currently supports Spheroid, Chebyshev, and Cylinder slices only")
    _linearization_variables_supported(variables(problem), canonical_variables) ||
        return LinearizationSupport(false,
                                    "EBCM analytical linearization supports unique canonical variables drawn from $(_ebcm_variable_list_message(canonical_variables))")

    return LinearizationSupport(true, "")
end

function _checked_ebcm_linearization_input(problem::LinearizationProblem,
                                          backend::EBCMLinearization,
                                          config)
    support = supports_linearization(problem, backend; output = :transition_matrix, config)
    Bool(support) ||
        throw(UnsupportedLinearization(backend, :transition_matrix, support.reason))
    return _ebcm_linearization_input(problem, config)
end

function linearize_transition_matrix(problem::LinearizationProblem,
                                     backend::EBCMLinearization; config = nothing)
    input = _checked_ebcm_linearization_input(problem, backend, config)
    𝐏s, 𝐔s = _ebcm_matrices(input)

    # 𝐐 = 𝐏 + i𝐔 per m-block depends only on 𝐏 and 𝐔, not on the differentiated
    # parameter. Factor each block once here and reuse it for the value block and
    # every Jacobian slice, instead of recomputing inv(𝐐) per parameter per block.
    factors = [_ebcm_factor(@. 𝐏 + 1im * 𝐔) for (𝐏, 𝐔) in zip(𝐏s, 𝐔s)]
    𝐓blocks = [-_ebcm_rdiv(𝐏, F) for (𝐏, F) in zip(𝐏s, factors)]
    value = _axisymmetric_transition_matrix_from_blocks(𝐓blocks)

    matrix_derivatives = _ebcm_matrix_derivatives(input, variables(problem))
    jacobian = map(matrix_derivatives) do (∂𝐏s, ∂𝐔s)
        ∂blocks = map(∂𝐏s, ∂𝐔s, 𝐓blocks, factors) do ∂𝐏, ∂𝐔, 𝐓, F
            ∂𝐐 = @. ∂𝐏 + 1im * ∂𝐔
            # ∂𝐓 = -(∂𝐏 + 𝐓 ∂𝐐) 𝐐⁻¹, reusing the block factor and value block.
            -_ebcm_rdiv(∂𝐏 + 𝐓 * ∂𝐐, F)
        end
        _axisymmetric_transition_matrix_from_blocks(∂blocks)
    end

    return LinearizationResult(value, jacobian, variables(problem);
                               metadata = (; backend = :ebcm_analytic,
                                           λ = input.λ,
                                           nₘₐₓ = input.nₘₐₓ,
                                           Ng = input.Ng))
end
