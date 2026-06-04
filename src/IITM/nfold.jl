@doc raw"""
```
transition_matrix_iitm(s::AbstractNFoldShape{N, T, CT}, λ, nₘₐₓ, Nr, Nϑ, Nφ; rₘᵢₙ) where {T, CT}
```

Use IITM to calculate the T-Matrix for a given scatterer and wavelength.

Parameters:

- `s`: the scatterer.
- `λ`: the wavelength.
- `nₘₐₓ`: the maximum order of the T-Matrix.
- `Nr`: the number of radial quadrature points to be used.
- `Nϑ`: the number of zenithal quadrature points to be used.
- `Nφ`: the number of azimuthal quadrature points to be used.

Keyword arguments:

- `rₘᵢₙ`: the starting point of the radial quadrature. Default to `rmin(s)`, which is the radius of the maximum inscribed sphere.

Returns:

- `𝐓`: an `AxisymmetricTransitionMatrix` struct representing the T-Matrix.
"""
function transition_matrix_iitm(s::AbstractNFoldShape{N, T, CT}, λ, nₘₐₓ, Nr, Nϑ,
        Nφ; rₘᵢₙ = rmin(s)) where {N, T, CT}
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

    if CT <: ComplexF64
        xφ = range(0, 2 * T(π) / N, length = Nφ + 1)[1:(end - 1)]
        wφ = 2 * T(π) / (N * Nφ)
    else
        xφ, wφ = gausslegendre(T, Nφ)
        @. xφ = (xφ + 1) * π / N
        @. wφ = π / N * wφ
    end

    fourier_workspace = CT <: ComplexF64 ? _azimuthal_fourier_workspace(Nφ, Nϑ) : nothing
    fourier_modes = CT <: ComplexF64 ? _azimuthal_fourier_mode_bins(nₘₐₓ, N) : nothing

    it = collect(OrderDegreeIterator(nₘₐₓ))
    its = [collect(enumerate([x for x in it if (x[2] % N + N) % N == i]))
           for i in 0:(N - 1)]
    𝐓s = [zeros(CT, 2 * length(its[i]), 2 * length(its[i])) for i in 1:N]

    # Initialize T-Matrix with Mie coefficients
    a, b = bhmie(T, k * rₘᵢₙ, s.m; nₘₐₓ = nₘₐₓ)
    for (it, 𝐓) in zip(its, 𝐓s)
        for (j, (n, _)) in it
            𝐓[2j - 1, 2j - 1] = -b[n]
            𝐓[2j, 2j] = -a[n]
        end
    end

    # Precalculate d, 𝜋 and τ
    d = OffsetArray(zeros(T, Nϑ, nₘₐₓ + 1, 2nₘₐₓ + 1), 1:Nϑ, 0:nₘₐₓ, (-nₘₐₓ):nₘₐₓ)
    𝜋 = similar(d)
    τ = similar(d)

    Threads.@threads for (i, m) in collect(Iterators.product(1:Nϑ, (-nₘₐₓ):nₘₐₓ))
        TransitionMatrices.wigner_d_recursion!(
            view(d, i, abs(m):nₘₐₓ, m), 0, m, nₘₐₓ, ϑ[i];
            deriv = view(τ, i, abs(m):nₘₐₓ, m))

        for n in max(abs(m), 1):nₘₐₓ
            𝜋[i, n, m] = TransitionMatrices.pi_func(T, m, n, ϑ[i]; d = d[i, n, m])
        end
    end

    # Precalculate coefficients
    a½ = [√(T(n * (n + 1))) for n in 1:nₘₐₓ]
    A = [√(T(2n + 1) / (2n * (n + 1))) for n in 1:nₘₐₓ]

    𝐉s = [zeros(CT, 3 * length(it), 2 * length(it)) for it in its]
    𝐇s = [zeros(CT, 3 * length(it), 2 * length(it)) for it in its]
    𝐆s = [zeros(CT, 3 * length(it), 3 * length(it)) for it in its]
    𝐔s = [zeros(CT, 3 * length(it), 3 * length(it)) for it in its]

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
        nₑₓₜᵣₐ = TransitionMatrices.estimate_ricattibesselj_extra_terms(nₘₐₓ, kr)
        TransitionMatrices.ricattibesselj!(ψ, ψ′, z, nₘₐₓ, nₑₓₜᵣₐ, kr)
        TransitionMatrices.ricattibessely!(χ, χ′, nₘₐₓ, kr)

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
        for (𝐉, 𝐇, 𝐆, it) in zip(𝐉s, 𝐇s, 𝐆s, its)
            for (i, (n, _)) in it
                𝐉[(3i - 2):(3i), (2i - 1):(2i)] .= 𝐉ᵈ[n]
                𝐇[(3i - 2):(3i), (2i - 1):(2i)] .= 𝐇ᵈ[n]

                # 𝐆 is averaged from both sides
                𝐆[(3i - 2):(3i), (3i - 2):(3i)] .= (𝐇ᵈ[n] * transpose(𝐉ᵈ[n]) .+
                                                    𝐉ᵈ[n] * transpose(𝐇ᵈ[n])) .*
                                                   (im * k / 2)
            end
        end

        # Calculate for each point whether it is within the scatterer
        ε = [refractive_index(s,
                 (r * sin(ϑ[i]) * cos(φ), r * sin(ϑ[i]) * sin(φ),
                     r * x[i]))^2 for φ in xφ, i in 1:Nϑ]
        fourier_coeffs = CT <: ComplexF64 ?
                         _azimuthal_fourier_coefficients(ε, nₘₐₓ, wφ,
            fourier_workspace,
            fourier_modes) : nothing

        # Calculate for each rem
        for (𝐓, 𝐉, 𝐇, 𝐆, 𝐔, it) in zip(𝐓s, 𝐉s, 𝐇s, 𝐆s, 𝐔s, its)
            Threads.@threads for (q, (n′, m′)) in it
                for (p, (n, m)) in it
                    sig = iseven(m + m′) ? 1 : -1

                    if has_symmetric_plane(s)
                        c = iseven(n + m + n′ + m′) ? 2 : 0
                        c̃ = 2 - c
                    else
                        c = 1
                        c̃ = 1
                    end

                    U = zero(SMatrix{3, 3, CT})

                    for i in 1:Nϑ
                        pptt = 𝜋[i, n, m] * 𝜋[i, n′, m′] + τ[i, n, m] * τ[i, n′, m′]
                        pttp = 𝜋[i, n, m] * τ[i, n′, m′] + τ[i, n, m] * 𝜋[i, n′, m′]
                        dd = d[i, n, m] * d[i, n′, m′]

                        if isnothing(fourier_coeffs)
                            for (j, φ) in enumerate(xφ)
                                ΔU = @SMatrix [c*pptt -c̃*im*pttp 0
                                               c̃*im*pttp c*pptt 0
                                               0 0 c * a½[n] * a½[n′] * dd/ε[j, i]]
                                U += w[i] * wφ[j] * cis((m′ - m) * φ) * (ε[j, i] - 1) * ΔU
                            end
                        else
                            coeff_ε, coeff_εinv = fourier_coeffs
                            freq = m′ - m
                            cε = coeff_ε[freq, i]
                            cεinv = coeff_εinv[freq, i]
                            U += w[i] * @SMatrix [c*pptt*cε -c̃*im*pttp*cε 0
                                           c̃*im*pttp*cε c*pptt*cε 0
                                           0 0 c * a½[n] * a½[n′] * dd * cεinv]
                        end
                    end

                    U *= (kr)^2 * sig * A[n] * A[n′] * N / 2 / π
                    𝐔[(3p - 2):(3p), (3q - 2):(3q)] .= U
                end
            end

            𝐐 = wri * _iitm_ldiv(𝐈 - wri * 𝐔 * 𝐆, 𝐔)
            𝐐ⱼⱼ = im * k * transpose(𝐉) * 𝐐 * 𝐉
            𝐐ⱼₕ = im * k * transpose(𝐉) * 𝐐 * 𝐇
            𝐐ₕⱼ = im * k * transpose(𝐇) * 𝐐 * 𝐉
            𝐐ₕₕ = im * k * transpose(𝐇) * 𝐐 * 𝐇

            # Eq (97) in Johnson (1988) Note: 𝐐ⱼₕ'==𝐐ₕⱼ only holds for spheres
            # Eq (38) in Bi et al. (2013)
            # Eq (2.40) in Doicu & Wriedt (2018)
            # Eq (5.71) in Hu (2018)
            # Eq (4.2.36) in Sun et al. (2019) Note: incorrect multiplication order
            𝐓 .= 𝐐ⱼⱼ + (𝐈 + 𝐐ⱼₕ) * _iitm_ldiv(𝐈 - 𝐓 * 𝐐ₕₕ, 𝐓) * (𝐈 + 𝐐ₕⱼ)
        end
    end

    𝐓′ = OffsetArray(zeros(CT, 2nₘₐₓ + 1, nₘₐₓ, 2nₘₐₓ + 1, nₘₐₓ, 2, 2), (-nₘₐₓ):nₘₐₓ,
        1:nₘₐₓ, (-nₘₐₓ):nₘₐₓ, 1:nₘₐₓ, 1:2, 1:2)

    for (𝐓, it) in zip(𝐓s, its)
        Threads.@threads for (j, (n′, m′)) in it
            for (i, (n, m)) in it
                𝐓′[m, n, m′, n′, 1, 1] = 𝐓[2i - 1, 2j - 1]
                𝐓′[m, n, m′, n′, 1, 2] = 𝐓[2i - 1, 2j]
                𝐓′[m, n, m′, n′, 2, 1] = 𝐓[2i, 2j - 1]
                𝐓′[m, n, m′, n′, 2, 2] = 𝐓[2i, 2j]
            end
        end
    end

    return TransitionMatrix{CT, nₘₐₓ, typeof(𝐓′)}(𝐓′)
end

@testitem "Spheroids can be solved by the nfold-shape solver" begin
    struct FourFoldSpheroid{T, CT} <: AbstractNFoldShape{4, T, CT}
        s::Spheroid{T, CT}
        m::CT
    end

    TransitionMatrices.rmin(s::FourFoldSpheroid) = rmin(s.s)
    TransitionMatrices.rmax(s::FourFoldSpheroid) = rmax(s.s)
    TransitionMatrices.refractive_index(s::FourFoldSpheroid, x) = refractive_index(s.s, x)

    s = Spheroid(1.0, 2.0, complex(1.5))
    ss = FourFoldSpheroid(s, s.m)

    nₘₐₓ = 5
    Nϑ = 50
    𝐓ₑ = calc_T(s, 2π, nₘₐₓ, Nϑ)

    Nr = 50
    Nφ = 100
    𝐓ᵢ = calc_T_iitm(ss, 2π, nₘₐₓ, Nr, 4Nϑ, Nφ)

    Cˢᶜᵃₑ = calc_Csca(𝐓ₑ)
    Cˢᶜᵃᵢ = calc_Csca(𝐓ᵢ)
    Cᵉˣᵗₑ = calc_Cext(𝐓ₑ)
    Cᵉˣᵗᵢ = calc_Cext(𝐓ᵢ)

    @test abs(Cˢᶜᵃₑ - Cˢᶜᵃᵢ) < 1e-2
    @test abs(Cᵉˣᵗₑ - Cᵉˣᵗᵢ) < 1e-2
end
