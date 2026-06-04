@doc raw"""
```
transition_matrix_iitm(s::AbstractAxisymmetricShape{T, CT}, λ, nₘₐₓ, Nr, Nϑ; rₘᵢₙ) where {T, CT}
```

Use IITM to calculate the T-Matrix for a given scatterer and wavelength.

Parameters:

- `s`: the axisymmetric scatterer.
- `λ`: the wavelength.
- `nₘₐₓ`: the maximum order of the T-Matrix.
- `Nr`: the number of radial quadrature points to be used.
- `Nϑ`: the number of zenithal quadrature points to be used.

Keyword arguments:

- `rₘᵢₙ`: the starting point of the radial quadrature. Default to `rmin(s)`, which is the radius of the maximum inscribed sphere.

Returns:

- `𝐓`: an `AxisymmetricTransitionMatrix` struct representing the T-Matrix.
"""
function transition_matrix_iitm(s::AbstractAxisymmetricShape{T, CT}, λ, nₘₐₓ, Nr, Nϑ
                                ; rₘᵢₙ = rmin(s)) where {T, CT}
    k = 2 * T(π) / λ
    rₘₐₓ = rmax(s)

    # Radial quadrature nodes and weights
    xr, wr = gausslegendre(T, Nr)
    @. xr = (rₘₐₓ - rₘᵢₙ) * (xr + 1) / 2 + rₘᵢₙ
    @. wr = (rₘₐₓ - rₘᵢₙ) / 2 * wr

    # Zenithal quadrature nodes and weights
    x, w = gausslegendre(T, Nϑ)
    ϑ = acos.(x)
    Nϑ = has_symmetric_plane(s) ? Nϑ ÷ 2 : Nϑ

    # Initialize T-Matrix with Mie coefficients
    a, b = bhmie(T, k * rₘᵢₙ, s.m; nₘₐₓ = nₘₐₓ)
    Ts = Matrix{CT}[]
    for m in 0:nₘₐₓ
        nn = nₘₐₓ + 1 - max(m, 1)
        offset = max(m, 1) - 1
        𝐓ₘ = zeros(CT, 2nn, 2nn)

        # Note that during calculation, we are using a different structure for 𝐓 and 𝐐
        for n in 1:nn
            𝐓ₘ[2n - 1, 2n - 1] = -b[n + offset]
            𝐓ₘ[2n, 2n] = -a[n + offset]
        end
        push!(Ts, 𝐓ₘ)
    end

    # Precalculate d, 𝜋 and τ
    d = OffsetArray(zeros(T, Nϑ, nₘₐₓ + 1, nₘₐₓ + 1), 1:Nϑ, 0:nₘₐₓ, 0:nₘₐₓ)
    𝜋 = similar(d)
    τ = similar(d)

    Threads.@threads for (i, m) in collect(Iterators.product(1:Nϑ, 0:nₘₐₓ))
        wigner_d_recursion!(view(d, i, m:nₘₐₓ, m), 0, m, nₘₐₓ, ϑ[i];
                            deriv = view(τ, i, m:nₘₐₓ, m))

        for n in max(m, 1):nₘₐₓ
            𝜋[i, n, m] = pi_func(T, m, n, ϑ[i]; d = d[i, n, m])
        end
    end

    # Precalculate coefficients
    a½ = [√(T(n * (n + 1))) for n in 1:nₘₐₓ]
    A = [√(T(2n + 1) / (2n * (n + 1))) for n in 1:nₘₐₓ]

    𝐉 = zeros(CT, 3nₘₐₓ, 2nₘₐₓ)
    𝐇 = zeros(CT, 3nₘₐₓ, 2nₘₐₓ)
    𝐆 = zeros(CT, 3nₘₐₓ, 3nₘₐₓ)

    nₑₓₜᵣₐ = TransitionMatrices.estimate_ricattibesselj_extra_terms(nₘₐₓ, k * rₘₐₓ)
    ψ = zeros(T, nₘₐₓ)
    z = zeros(T, nₘₐₓ + nₑₓₜᵣₐ)
    ψ′ = similar(ψ)
    χ = similar(ψ)
    χ′ = similar(ψ)

    # Radial recursion
    for (r, wri) in zip(xr, wr)
        @debug "Calculating layer r = $r"

        # Calculate Ricatti-Bessel functions and derivatives
        kr = k * r
        nₑₓₜᵣₐ = estimate_ricattibesselj_extra_terms(nₘₐₓ, kr)
        ricattibesselj!(ψ, ψ′, z, nₘₐₓ, nₑₓₜᵣₐ, kr)
        ricattibessely!(χ, χ′, nₘₐₓ, kr)

        # Since we use Ricatti-Bessel instead of spherical Bessel functions,
        # we need to divide the values by an extra `kr`
        𝐉ᵈ = [@SMatrix [ψ[n]/kr 0
                        0 ψ′[n]/kr
                        0 a½[n] * ψ[n]/(kr)^2] for n in 1:nₘₐₓ]
        𝐘ᵈ = [@SMatrix [χ[n]/kr 0
                        0 χ′[n]/kr
                        0 a½[n] * χ[n]/(kr)^2] for n in 1:nₘₐₓ]
        𝐇ᵈ = @. 𝐉ᵈ + 1im * 𝐘ᵈ

        # Make block diagonal matrices
        for n in 1:nₘₐₓ
            𝐉[(3n - 2):(3n), (2n - 1):(2n)] .= 𝐉ᵈ[n]
            𝐇[(3n - 2):(3n), (2n - 1):(2n)] .= 𝐇ᵈ[n]

            # 𝐆 is averaged from both sides
            𝐆[(3n - 2):(3n), (3n - 2):(3n)] .= (𝐇ᵈ[n] * transpose(𝐉ᵈ[n]) .+
                                                𝐉ᵈ[n] * transpose(𝐇ᵈ[n])) .* (im * k / 2)
        end

        # Calculate for each point whether it is within the scatterer
        within = [(r * sin(ϑ[i]), 0, r * x[i]) ∈ s for i in 1:Nϑ]
        ε = [within[i] ? s.m^2 : one(CT) for i in eachindex(within)]

        # Calculate for each m
        Threads.@threads for m in 0:nₘₐₓ
            nₘᵢₙ = max(1, m)
            nn = nₘₐₓ - nₘᵢₙ + 1
            𝐔 = zeros(CT, 3nn, 3nn)

            for n in nₘᵢₙ:nₘₐₓ, n′ in nₘᵢₙ:nₘₐₓ
                U = zero(SMatrix{3, 3, CT})
                if has_symmetric_plane(s)
                    c = iseven(n + n′) ? 2 : 0
                    c̃ = 2 - c
                else
                    c = 1
                    c̃ = 1
                end

                for i in 1:Nϑ
                    pptt = 𝜋[i, n, m] * 𝜋[i, n′, m] + τ[i, n, m] * τ[i, n′, m]
                    pttp = 𝜋[i, n, m] * τ[i, n′, m] + τ[i, n, m] * 𝜋[i, n′, m]
                    dd = d[i, n, m] * d[i, n′, m]
                    ΔU = @SMatrix [c*pptt -c̃*im*pttp 0
                                   c̃*im*pttp c*pptt 0
                                   0 0 c * a½[n] * a½[n′] * dd/ε[i]]
                    U += w[i] * (ε[i] - 1) * ΔU
                end

                U *= A[n] * A[n′] * (kr)^2
                view(𝐔, (3(n - nₘᵢₙ + 1) - 2):(3(n - nₘᵢₙ + 1)),
                (3(n′ - nₘᵢₙ + 1) - 2):(3(n′ - nₘᵢₙ + 1))) .= U
            end

            𝐉ᵥ = view(𝐉, (3nₘᵢₙ - 2):(3nₘₐₓ), (2nₘᵢₙ - 1):(2nₘₐₓ))
            𝐇ᵥ = view(𝐇, (3nₘᵢₙ - 2):(3nₘₐₓ), (2nₘᵢₙ - 1):(2nₘₐₓ))
            𝐆ᵥ = view(𝐆, (3nₘᵢₙ - 2):(3nₘₐₓ), (3nₘᵢₙ - 2):(3nₘₐₓ))
            𝐐 = wri * _iitm_ldiv(𝐈 - wri * 𝐔 * 𝐆ᵥ, 𝐔)

            𝐐ⱼⱼ = im * k * transpose(𝐉ᵥ) * 𝐐 * 𝐉ᵥ
            𝐐ⱼₕ = im * k * transpose(𝐉ᵥ) * 𝐐 * 𝐇ᵥ
            𝐐ₕⱼ = im * k * transpose(𝐇ᵥ) * 𝐐 * 𝐉ᵥ
            𝐐ₕₕ = im * k * transpose(𝐇ᵥ) * 𝐐 * 𝐇ᵥ

            # Eq (97) in Johnson (1988) Note: 𝐐ⱼₕ'==𝐐ₕⱼ only holds for spheres
            # Eq (38) in Bi et al. (2013)
            # Eq (2.40) in Doicu & Wriedt (2018)
            # Eq (5.71) in Hu (2018)
            # Eq (4.2.36) in Sun et al. (2019) Note: incorrect multiplication order
            Ts[m + 1] = 𝐐ⱼⱼ +
                        (𝐈 + 𝐐ⱼₕ) * _iitm_ldiv(𝐈 - Ts[m + 1] * 𝐐ₕₕ, Ts[m + 1]) *
                        (𝐈 + 𝐐ₕⱼ)
        end
    end

    Threads.@threads for m in 0:nₘₐₓ
        T′ = similar(Ts[m + 1])
        nₘᵢₙ = max(1, m)
        nn = nₘₐₓ - nₘᵢₙ + 1

        for n in 1:nn
            for n′ in 1:nn
                T′[n, n′] = Ts[m + 1][2n - 1, 2n′ - 1]
                T′[n, nn + n′] = Ts[m + 1][2n - 1, 2n′]
                T′[nn + n, n′] = Ts[m + 1][2n, 2n′ - 1]
                T′[nn + n, nn + n′] = Ts[m + 1][2n, 2n′]
            end
        end

        Ts[m + 1] = T′
    end

    return AxisymmetricTransitionMatrix{CT, nₘₐₓ, typeof(Ts), T}(Ts)
end

@testitem "IITM is correct for spheres" begin
    using TransitionMatrices: transition_matrix_iitm

    r = 5.0
    λ = 2π
    k = 2π / λ
    m = 1.311 + 0.01im
    s = Spheroid(r, r, m)

    nₘₐₓ = 10
    𝐓ₘ = MieTransitionMatrix{ComplexF64, nₘₐₓ}(k * r, m)

    Nr = 50
    Nϑ = 100
    𝐓ᵢ = transition_matrix_iitm(s, λ, nₘₐₓ, Nr, Nϑ; rₘᵢₙ = 4.0)

    Cˢᶜᵃₘ = calc_Csca(𝐓ₘ)
    Cˢᶜᵃᵢ = calc_Csca(𝐓ᵢ)
    Cᵉˣᵗₘ = calc_Cext(𝐓ₘ)
    Cᵉˣᵗᵢ = calc_Cext(𝐓ᵢ)

    @test abs(Cˢᶜᵃₘ - Cˢᶜᵃᵢ) < 1e-2
    @test abs(Cᵉˣᵗₘ - Cᵉˣᵗᵢ) < 1e-2
end

@testitem "IITM is correct for spheroids" begin
    using TransitionMatrices: transition_matrix_iitm

    a = 5.0
    c = 6.0
    λ = 2π
    k = 2π / λ
    m = 1.311 + 0.01im
    s = Spheroid(a, c, m)
    nₘₐₓ = 15
    Ng = 40
    𝐓ₑ = transition_matrix(s, λ, nₘₐₓ, Ng)

    Nr = 500
    𝐓ᵢ = transition_matrix_iitm(s, λ, nₘₐₓ, Nr, 4Ng)

    Cˢᶜᵃₑ = calc_Csca(𝐓ₑ)
    Cˢᶜᵃᵢ = calc_Csca(𝐓ᵢ)
    Cᵉˣᵗₑ = calc_Cext(𝐓ₑ)
    Cᵉˣᵗᵢ = calc_Cext(𝐓ᵢ)

    @test abs(Cˢᶜᵃₑ - Cˢᶜᵃᵢ) < 1e-2
    @test abs(Cᵉˣᵗₑ - Cᵉˣᵗᵢ) < 1e-2
end
