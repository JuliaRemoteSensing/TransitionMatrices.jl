function estimate_ricattibesselj_extra_terms(n, x::Real)
    tb = max(n, abs(x))
    ceil(Int, tb - n + 8.0 * √tb + 3)
end

function estimate_ricattibesselj_extra_terms(n, x)
    tb = max(n, abs(x))
    ceil(Int, tb - n + 4.0 * ∛tb + 8.0 * √tb + 5)
end

@doc raw"""
```
ricattibesselj!(ψₙ, ψₙ′, z, nₘₐₓ, nₑₓₜᵣₐ, x)
```

Parameters:

- `ψ`: output array for ``\psi_n(x)``.
- `ψ′`: output array for ``\psi^{\prime}_n(x)``.
- `z`: auxiliary array for downward recursion.
- `nₘₐₓ`: maximum order of ``\psi_n(x)``.
- `nₑₓₜᵣₐ`: extra terms for downward recursion of ``z``.
- `x`: argument.

In-place version of [`ricattibesselj`](@ref).
"""
function ricattibesselj!(ψ, ψ′, z, nₘₐₓ, nₑₓₜᵣₐ, x)
    nₜₒₜₐₗ = nₘₐₓ + nₑₓₜᵣₐ
    x⁻¹ = one(x) / x
    if length(z) < nₜₒₜₐₗ
        resize!(z, nₜₒₜₐₗ)
    end
    z[nₜₒₜₐₗ] = x / (2nₜₒₜₐₗ + 1)
    for n in (nₜₒₜₐₗ - 1):-1:1
        z[n] = one(x) / ((2n + 1) * x⁻¹ - z[n + 1])
    end
    z₀ = one(x) / (x⁻¹ - z[1])
    y₀ = z₀ * cos(x)
    ψ[1] = y₀ * z[1]
    ψ′[1] = y₀ - ψ[1] * x⁻¹
    for n in 2:nₘₐₓ
        ψ[n] = ψ[n - 1] * z[n]
        ψ′[n] = ψ[n - 1] - ψ[n] * x⁻¹ * n
    end
end

@doc raw"""
```
ricattibesselj(nₘₐₓ::Integer, nₑₓₜᵣₐ::Integer, x)
```

Calculate Ricatti-Bessel function of the first kind, ``\psi_n(x)``, and its derivative ``\psi^{\prime}_n(x)`` for ``1\le n\le n_{\max}``.

``\psi_n(x)`` is defined as

```math
\psi_n(x) = xj_n(x)
```

where ``j_n(x)`` is the Bessel function of the first kind.

Parameters:

- `nₘₐₓ`: maximum order of ``\psi_n(x)``.
- `nₑₓₜᵣₐ`: extra terms for downward recursion of ``z``.
- `x`: argument.

Returns: (`ψ`, `ψ′`)

- `ψ`: array of ``\psi_n(x)``.
- `ψ′`: array of ``\psi^{\prime}_n(x)``.
"""
function ricattibesselj(nₘₐₓ::Integer, nₑₓₜᵣₐ::Integer, x)
    T = typeof(x)
    ψ = zeros(T, nₘₐₓ)
    ψ′ = zeros(T, nₘₐₓ)
    z = zeros(T, nₘₐₓ + nₑₓₜᵣₐ)
    ricattibesselj!(ψ, ψ′, z, nₘₐₓ, nₑₓₜᵣₐ, x)
    return ψ, ψ′
end

@testitem "Ricatti-Bessel ψ and ψ′" begin
    # Check correctness of real-value Ricatti-Bessel functions against `Bessels.jl`.

    using TransitionMatrices: ricattibesselj, estimate_ricattibesselj_extra_terms
    using GSL: sf_bessel_jl_array

    @testset "x = $x, n ∈ [1, $n]" for x in (0.01, 0.1, 1.0, 10.0, 100.0, 1000.0),
                                       n in (40,)

        ψ, ψ′ = ricattibesselj(n, estimate_ricattibesselj_extra_terms(n, x), x)

        j = sf_bessel_jl_array(n, x)
        j′ = [j[i] - (i + 1) / x * j[i + 1] for i in 1:n]
        ψ₁ = j[2:end] .* x
        ψ₁′ = j′ .* x .+ j[2:end]

        @test all(ψ .≈ ψ₁)
        @test all(ψ′ .≈ ψ₁′)
    end
end

@doc raw"""
```
ricattibessely!(χ, χ′, nₘₐₓ, x)
```

Parameters:

- `χ`: output array for ``\chi_n(x)``.
- `χ′`: output array for ``\chi^{\prime}_n(x)``.
- `nₘₐₓ`: maximum order of ``\chi_n(x)``.
- `x`: argument.

In-place version of [`ricattibessely`](@ref).
"""
function ricattibessely!(χ, χ′, nₘₐₓ, x)
    x⁻¹ = one(x) / x
    χ[1] = -cos(x) * x⁻¹ - sin(x)
    χ[2] = (-3x⁻¹^2 + 1) * cos(x) - 3x⁻¹ * sin(x)
    for n in 2:(nₘₐₓ - 1)
        χ[n + 1] = (2n + 1) * x⁻¹ * χ[n] - χ[n - 1]
    end
    χ′[1] = -(cos(x) + x⁻¹ * χ[1])
    for n in 2:nₘₐₓ
        χ′[n] = χ[n - 1] - n * x⁻¹ * χ[n]
    end
end

@doc raw"""
```
ricattibessely(nₘₐₓ::Integer, x)
```

Calculate Ricatti-Bessel function of the second kind, ``\chi_n(x)``, and its derivative ``\chi^{\prime}_n(x)`` for ``1\le n\le n_{\max}``.

``\chi_n(x)`` is defined as

```math
\chi_n(x) = xy_n(x)
```

where ``y_n(x)`` is the Bessel function of the second kind.

Parameters:

- `nₘₐₓ`: maximum order of ``\chi_n(x)``.
- `x`: argument.

Returns: (`χ`, `χ′`)

- `χ`: Array of ``\chi_n(x)``.
- `χ′`: Array of ``\chi^{\prime}_n(x)``.
"""
function ricattibessely(nₘₐₓ, x)
    T = typeof(x)
    χ = zeros(T, nₘₐₓ)
    χ′ = zeros(T, nₘₐₓ)
    ricattibessely!(χ, χ′, nₘₐₓ, x)
    return χ, χ′
end

@testitem "Ricatti-Bessel χ and χ′" begin
    # Check correctness of real-value Ricatti-Bessel functions against GSL.

    using TransitionMatrices: ricattibessely
    using GSL: sf_bessel_yl_array

    @testset "x = $x, n ∈ [1, $n]" for x in (0.1, 1.0, 10.0, 100.0), n in (40,)
        χ, χ′ = ricattibessely(n, x)

        y = sf_bessel_yl_array(n, x)
        y′ = [y[i] - (i + 1) / x * y[i + 1] for i in 1:n]
        χ₁ = y[2:end] .* x
        χ₁′ = y′ .* x .+ y[2:end]

        @test all(χ .≈ χ₁)
        @test all(χ′ .≈ χ₁′)
    end
end
