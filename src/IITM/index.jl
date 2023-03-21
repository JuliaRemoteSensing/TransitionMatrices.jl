Base.@kwdef struct IITMConfiguration
    N::Int
    Ng::Int
    Nr::Int
end

function transition_matrix(s::AbstractAxisymmetricShape{T, CT}, λ,
                           cfg::IITMConfiguration) where {T, CT}
    (; N, Ng, Nr) = cfg
    k = 2π / λ
    rₘᵢₙ = rmin(s) - 0.5
    rₘₐₓ = rmax(s)
    a, b = bhmie(T, k * rₘᵢₙ, s.m; nₘₐₓ = N)

    xr, wr = gausslegendre(Nr)
    @. xr = (rₘₐₓ - rₘᵢₙ) * (xr + 1) / 2 + rₘᵢₙ
    @. wr = (rₘₐₓ - rₘᵢₙ) / 2 * wr

    x, w = gausslegendre(Ng)
    ϑ = acos.(x)

    Ts = Matrix{CT}[]
    for m in 0:N
        n = N + 1 - max(m, 1)
        offset = max(m, 1) - 1
        𝐓ₘ = zeros(CT, 2n, 2n)

        # Note that during calculation we are using a different structure for 𝐓 and 𝐐
        for i in 1:n
            𝐓ₘ[2i - 1, 2i - 1] = -b[i + offset]
            𝐓ₘ[2i, 2i] = -a[i + offset]
        end
        push!(Ts, 𝐓ₘ)
    end

    d = OffsetArray(zeros(T, Ng, N + 1, N + 1), 1:Ng, 0:N, 0:N)
    𝜋 = similar(d)
    τ = similar(d)

    Threads.@threads for i in eachindex(ϑ)
        for m in 0:N
            wigner_d_recursion!(view(d, i, m:N, m), 0, m, N, ϑ[i];
                                deriv = view(τ, i, m:N, m))

            for n in max(m, 1):N
                𝜋[i, n, m] = pi_func(T, m, n, ϑ[i]; d = d[i, n, m])
            end
        end
    end

    as = [√(T(n * (n + 1))) for n in 1:N]
    A = [√(T(2n + 1) / (2n * (n + 1))) for n in 1:N]
    for (r, wri) in zip(xr, wr)
        @debug "Calculating layer r = $r"
        kr = k * r
        nₑₓₜᵣₐ = estimate_ricattibesselj_extra_terms(N, kr)
        ψ = zeros(T, N)
        z = zeros(T, N + nₑₓₜᵣₐ)
        ψ′ = similar(ψ)
        χ = similar(ψ)
        χ′ = similar(ψ)
        ricattibesselj!(ψ, ψ′, z, N, nₑₓₜᵣₐ, kr)
        ricattibessely!(χ, χ′, N, kr)

        𝐉 = [@SMatrix [ψ[n]/kr 0
                       0 ψ′[n]/kr
                       0 as[n] * ψ[n]/kr] for n in 1:N]
        𝐋 = [@SMatrix [χ[n]/kr 0
                       0 χ′[n]/kr
                       0 as[n] * χ[n]/kr] for n in 1:N]
        𝐇 = @. 𝐉 + 1im * 𝐋
        𝐆 = [(H * J' + J * H') * (im * k / 2) for (J, H) in zip(𝐉, 𝐇)]

        𝐉 = collect(mortar(Diagonal(𝐉)))
        𝐇 = collect(mortar(Diagonal(𝐇)))
        𝐆 = collect(mortar(Diagonal(𝐆)))

        within = [(r * sin(ϑ[i]), 0, r * x[i]) ∈ s for i in eachindex(ϑ)]
        𝑚² = [within[i] ? s.m^2 : one(CT) for i in eachindex(within)]

        for m in 0:N
            nₘᵢₙ = max(1, m)
            nn = N - nₘᵢₙ + 1
            𝐔 = zeros(CT, 3nn, 3nn)

            for (n, n′) in collect(Iterators.product(nₘᵢₙ:N, nₘᵢₙ:N))
                U = zero(SMatrix{3, 3, CT})
                for i in eachindex(ϑ)
                    pptt = 𝜋[i, n, m] * 𝜋[i, n′, m] + τ[i, n, m] * τ[i, n′, m]
                    pttp = 𝜋[i, n, m] * τ[i, n′, m] + τ[i, n, m] * 𝜋[i, n′, m]
                    dd = d[i, n, m] * d[i, n′, m]
                    ΔU = @SMatrix [pptt -im*pttp 0
                                   im*pttp pptt 0
                                   0 0 as[n] * as[n′] * dd/𝑚²[i]]
                    U += ΔU * w[i] * (𝑚²[i] - 1)
                end

                U *= (kr)^2 * A[n] * A[n′] / 2
                view(𝐔, (3(n - nₘᵢₙ + 1) - 2):(3(n - nₘᵢₙ + 1)),
                (3(n′ - nₘᵢₙ + 1) - 2):(3(n′ - nₘᵢₙ + 1))) .= U
            end

            𝐆ᵥ = view(𝐆, (3nₘᵢₙ - 2):(3N), (3nₘᵢₙ - 2):(3N))
            𝐐 = wri * inv(𝐈 - wri * 𝐔 * 𝐆ᵥ) * 𝐔

            𝐇ᵥ = view(𝐇, (3nₘᵢₙ - 2):(3N), (2nₘᵢₙ - 1):(2N))
            𝐉ᵥ = view(𝐉, (3nₘᵢₙ - 2):(3N), (2nₘᵢₙ - 1):(2N))
            𝐐ₕₕ = 𝐇ᵥ' * 𝐐 * 𝐇ᵥ
            𝐐ₕⱼ = 𝐇ᵥ' * 𝐐 * 𝐉ᵥ
            𝐐ⱼₕ = 𝐉ᵥ' * 𝐐 * 𝐇ᵥ
            𝐐ⱼⱼ = 𝐉ᵥ' * 𝐐 * 𝐉ᵥ

            Ts[m + 1] = 𝐐ⱼⱼ +
                        (𝐈 + 𝐐ⱼₕ) * inv(𝐈 - Ts[m + 1] * 𝐐ₕₕ) * Ts[m + 1] *
                        (𝐈 + 𝐐ₕⱼ)
        end
    end

    for m in 0:N
        T′ = similar(Ts[m + 1])
        nn = N - max(1, m) + 1

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

    return AxisymmetricTransitionMatrix{CT, N, typeof(Ts), T}(Ts)
end
