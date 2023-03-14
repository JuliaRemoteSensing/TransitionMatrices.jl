# === Some functions that are not used in the package, but are kept for reference ===

function ricattibesselj!(ψ::AbstractVector{<:ArbLike}, ψ′, z, nₘₐₓ, nₑₓₜᵣₐ, x)
    ψ₀ = Arb(0)
    half = Arb(1 // 2)
    coeff = √(2x * π) * x
    Arblib.hypgeom_bessel_j!(ψ₀, half, x)
    ψ₀ *= coeff
    for n in 1:nₘₐₓ
        Arblib.hypgeom_bessel_j!(ψ[n], half + n, x)
        ψ[n] *= coeff
    end

    x⁻¹ = one(x) / x
    ψ′[1] = ψ₀ - ψ[1] * x⁻¹
    for n in 2:nₘₐₓ
        ψ′[n] = ψ[n - 1] - ψ[n] * x⁻¹ * n
    end
end

function ricattibesselj!(ψ::AbstractVector{<:AcbLike}, ψ′, z, nₘₐₓ, nₑₓₜᵣₐ, x)
    ψ₀ = Acb(0)
    half = Acb(1 // 2)
    coeff = √(2x * π) * x
    Arblib.hypgeom_bessel_j!(ψ₀, half, x)
    ψ₀ *= coeff
    for n in 1:nₘₐₓ
        Arblib.hypgeom_bessel_j!(ψ[n], half + n, x)
        ψ[n] *= coeff
    end

    x⁻¹ = one(x) / x
    ψ′[1] = ψ₀ - ψ[1] * x⁻¹
    for n in 2:nₘₐₓ
        ψ′[n] = ψ[n - 1] - ψ[n] * x⁻¹ * n
    end
end

function ricattibessely!(χ::AbstractVector{<:ArbLike}, χ′, nₘₐₓ, x)
    χ₀ = Arb(0)
    half = Arb(1 // 2)
    coeff = √(2x * π) * x
    Arblib.hypgeom_bessel_y!(χ₀, half, x)
    χ₀ *= coeff
    for n in 1:nₘₐₓ
        Arblib.hypgeom_bessel_y!(χ[n], half + n, x)
        χ[n] *= coeff
    end

    x⁻¹ = one(x) / x
    χ′[1] = χ₀ - χ[1] * x⁻¹
    for n in 2:nₘₐₓ
        χ′[n] = χ[n - 1] - χ[n] * x⁻¹ * n
    end
end
