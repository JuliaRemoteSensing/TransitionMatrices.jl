const _IITM_FIXED_GEOMETRY_LINEARIZATION_VARIABLES = (:mᵣ, :mᵢ, :λ)

function _iitm_variable_list_message(canonical)
    return join((":" * String(var) for var in canonical), ", ", " and ")
end

function _iitm_linearization_input(problem::LinearizationProblem, config,
        x = problem.x)
    rebuilt = rebuild(problem, x)
    shape = _linearization_property(config, :shape;
        default = _linearization_property(rebuilt, :shape))
    λ = _linearization_property(config, :λ;
        default = _linearization_property(rebuilt, :λ))
    nₘₐₓ = _linearization_property(config, :nₘₐₓ;
        default = _linearization_property(rebuilt, :nₘₐₓ))
    Nr = _linearization_property(config, :Nr;
        default = _linearization_property(rebuilt, :Nr))
    Nϑ = _linearization_property(config, :Nϑ;
        default = _linearization_property(rebuilt, :Nϑ))
    Nφ = _linearization_property(config, :Nφ;
        default = _linearization_property(rebuilt, :Nφ))

    if isnothing(shape) || isnothing(λ) || isnothing(nₘₐₓ) ||
       isnothing(Nr) || isnothing(Nϑ)
        return nothing
    end

    return (; shape, λ, nₘₐₓ = Int(nₘₐₓ), Nr = Int(Nr), Nϑ = Int(Nϑ),
        Nφ = isnothing(Nφ) ? nothing : Int(Nφ))
end

function _iitm_real_type(shape, λ)
    return typeof(float(real(λ)) + float(real(shape.m)) + float(rmin(shape)) +
                  float(rmax(shape)))
end

function _iitm_complex_refractive_index(::Type{RT}, m) where {RT}
    return complex(RT(real(m)), RT(imag(m)))
end

function _iitm_fixed_geometry_parameter_derivatives(input, variables)
    RT = _iitm_real_type(input.shape, input.λ)
    CT = complex(RT)
    k = 2 * RT(π) / RT(input.λ)
    m = _iitm_complex_refractive_index(RT, input.shape.m)
    ∂k = zeros(RT, length(variables))
    ∂m = zeros(CT, length(variables))

    for (j, variable) in enumerate(variables)
        if variable == :λ
            ∂k[j] = -k / RT(input.λ)
        elseif variable == :mᵣ
            ∂m[j] = one(CT)
        elseif variable == :mᵢ
            ∂m[j] = one(RT) * im
        end
    end

    return RT, CT, k, m, ∂k, ∂m
end

_iitm_nonzero_derivative_indices(∂x) = [j for j in eachindex(∂x) if !iszero(∂x[j])]

function _iitm_ricatti_argument_derivatives!(∂f, ∂f′, f, f′, z, ∂z)
    for n in eachindex(f, f′)
        ∂f[n] = f′[n] * ∂z
        ∂f′[n] = (n * (n + 1) / z^2 - 1) * f[n] * ∂z
    end
end

function _iitm_fill_radial_blocks!(𝐉, 𝐇, 𝐆, ∂𝐉s, ∂𝐇s, ∂𝐆s,
        ψ, ψ′, χ, χ′, ∂ψs, ∂ψ′s, ∂χs, ∂χ′s,
        kr, ∂kr, k, ∂k, a½, nₘₐₓ,
        wavenumber_indices = eachindex(∂k))
    CT = eltype(𝐉)
    fill!(𝐉, zero(CT))
    fill!(𝐇, zero(CT))
    fill!(𝐆, zero(CT))
    for ∂𝐉 in ∂𝐉s
        fill!(∂𝐉, zero(CT))
    end
    for ∂𝐇 in ∂𝐇s
        fill!(∂𝐇, zero(CT))
    end
    for ∂𝐆 in ∂𝐆s
        fill!(∂𝐆, zero(CT))
    end

    for n in 1:nₘₐₓ
        𝐉ᵈ = @SMatrix [ψ[n]/kr 0
                       0 ψ′[n]/kr
                       0 a½[n] * ψ[n]/kr^2]
        𝐘ᵈ = @SMatrix [χ[n]/kr 0
                       0 χ′[n]/kr
                       0 a½[n] * χ[n]/kr^2]
        𝐇ᵈ = 𝐉ᵈ + 1im * 𝐘ᵈ

        𝐉[(3n - 2):(3n), (2n - 1):(2n)] .= 𝐉ᵈ
        𝐇[(3n - 2):(3n), (2n - 1):(2n)] .= 𝐇ᵈ

        𝐆ᵈ_base = 𝐇ᵈ * transpose(𝐉ᵈ) + 𝐉ᵈ * transpose(𝐇ᵈ)
        𝐆[(3n - 2):(3n), (3n - 2):(3n)] .= 𝐆ᵈ_base .* (im * k / 2)

        for j in wavenumber_indices
            ∂ψ = ∂ψs[j][n]
            ∂ψ′ = ∂ψ′s[j][n]
            ∂χ = ∂χs[j][n]
            ∂χ′ = ∂χ′s[j][n]
            ∂krⱼ = ∂kr[j]

            ∂𝐉ᵈ = @SMatrix [∂ψ/kr - ψ[n] * ∂krⱼ/kr^2 0
                            0 ∂ψ′/kr - ψ′[n] * ∂krⱼ/kr^2
                            0 a½[n] * (∂ψ/kr^2 - 2ψ[n] * ∂krⱼ/kr^3)]
            ∂𝐘ᵈ = @SMatrix [∂χ/kr - χ[n] * ∂krⱼ/kr^2 0
                            0 ∂χ′/kr - χ′[n] * ∂krⱼ/kr^2
                            0 a½[n] * (∂χ/kr^2 - 2χ[n] * ∂krⱼ/kr^3)]
            ∂𝐇ᵈ = ∂𝐉ᵈ + 1im * ∂𝐘ᵈ
            ∂𝐆ᵈ_base = ∂𝐇ᵈ * transpose(𝐉ᵈ) + 𝐇ᵈ * transpose(∂𝐉ᵈ) +
                       ∂𝐉ᵈ * transpose(𝐇ᵈ) + 𝐉ᵈ * transpose(∂𝐇ᵈ)
            ∂𝐆ᵈ = ∂𝐆ᵈ_base .* (im * k / 2) .+
                  𝐆ᵈ_base .* (im * ∂k[j] / 2)

            ∂𝐉s[j][(3n - 2):(3n), (2n - 1):(2n)] .= ∂𝐉ᵈ
            ∂𝐇s[j][(3n - 2):(3n), (2n - 1):(2n)] .= ∂𝐇ᵈ
            ∂𝐆s[j][(3n - 2):(3n), (3n - 2):(3n)] .= ∂𝐆ᵈ
        end
    end

    return nothing
end

function _iitm_axisymmetric_shell_u!(𝐔, ∂𝐔s, s, r, x, w, ε, within, d, 𝜋,
        τ, a½, A, k, ∂k, m, ∂m, nₘₐₓ, m_order,
        ∂Us = nothing,
        material_indices = _iitm_nonzero_derivative_indices(∂m))
    nₘᵢₙ = max(1, m_order)
    p = length(∂m)
    CT = typeof(m)
    fill!(𝐔, zero(CT))
    for ∂𝐔 in ∂𝐔s
        fill!(∂𝐔, zero(CT))
    end
    kr = k * r
    if isnothing(∂Us)
        ∂Us = [zero(SMatrix{3, 3, CT}) for _ in 1:p]
    end

    for n in nₘᵢₙ:nₘₐₓ, n′ in nₘᵢₙ:nₘₐₓ

        U = zero(SMatrix{3, 3, CT})
        fill!(∂Us, zero(SMatrix{3, 3, CT}))
        if has_symmetric_plane(s)
            c = iseven(n + n′) ? 2 : 0
            c̃ = 2 - c
        else
            c = 1
            c̃ = 1
        end

        for i in eachindex(x)
            pptt = 𝜋[i, n, m_order] * 𝜋[i, n′, m_order] +
                   τ[i, n, m_order] * τ[i, n′, m_order]
            pttp = 𝜋[i, n, m_order] * τ[i, n′, m_order] +
                   τ[i, n, m_order] * 𝜋[i, n′, m_order]
            dd = d[i, n, m_order] * d[i, n′, m_order]
            ΔU = @SMatrix [c*pptt -c̃*im*pttp 0
                           c̃*im*pttp c*pptt 0
                           0 0 c * a½[n] * a½[n′] * dd/ε[i]]
            contrast = ε[i] - 1
            U += w[i] * contrast * ΔU

            if within[i]
                for j in material_indices
                    ∂ε = 2m * ∂m[j]
                    ∂ΔU = @SMatrix [zero(CT) zero(CT) zero(CT)
                                    zero(CT) zero(CT) zero(CT)
                                    zero(CT) zero(CT) -c * a½[n] * a½[n′] * dd * ∂ε/ε[i]^2]
                    ∂Us[j] += w[i] * (∂ε * ΔU + contrast * ∂ΔU)
                end
            end
        end

        prefactor = A[n] * A[n′] * kr^2
        U_unscaled = U
        U = prefactor * U_unscaled
        for j in 1:p
            ∂prefactor = A[n] * A[n′] * 2kr * ∂k[j] * r
            ∂Us[j] = ∂prefactor * U_unscaled + prefactor * ∂Us[j]
        end

        rows = (3(n - nₘᵢₙ + 1) - 2):(3(n - nₘᵢₙ + 1))
        cols = (3(n′ - nₘᵢₙ + 1) - 2):(3(n′ - nₘᵢₙ + 1))
        𝐔[rows, cols] .= U
        for j in 1:p
            ∂𝐔s[j][rows, cols] .= ∂Us[j]
        end
    end

    return 𝐔, ∂𝐔s
end

function _iitm_axisymmetric_shell_u(s, r, x, w, ε, within, d, 𝜋, τ, a½, A,
        k, ∂k, m, ∂m, nₘₐₓ, m_order)
    nₘᵢₙ = max(1, m_order)
    nn = nₘₐₓ - nₘᵢₙ + 1
    CT = typeof(m)
    𝐔 = zeros(CT, 3nn, 3nn)
    ∂𝐔s = [zeros(CT, 3nn, 3nn) for _ in eachindex(∂m)]
    return _iitm_axisymmetric_shell_u!(𝐔, ∂𝐔s, s, r, x, w, ε, within, d,
        𝜋, τ, a½, A, k, ∂k, m, ∂m, nₘₐₓ,
        m_order)
end

function _iitm_project_q_block_values(k, 𝐉, 𝐇, 𝐐)
    𝐉ᵀ𝐐 = transpose(𝐉) * 𝐐
    𝐇ᵀ𝐐 = transpose(𝐇) * 𝐐
    𝐐ⱼⱼ = im * k * 𝐉ᵀ𝐐 * 𝐉
    𝐐ⱼₕ = im * k * 𝐉ᵀ𝐐 * 𝐇
    𝐐ₕⱼ = im * k * 𝐇ᵀ𝐐 * 𝐉
    𝐐ₕₕ = im * k * 𝐇ᵀ𝐐 * 𝐇

    return 𝐐ⱼⱼ, 𝐐ⱼₕ, 𝐐ₕⱼ, 𝐐ₕₕ, 𝐉ᵀ𝐐, 𝐇ᵀ𝐐
end

function _iitm_project_q_blocks(k, ∂k, 𝐉, 𝐇, 𝐐, ∂𝐉, ∂𝐇, ∂𝐐)
    𝐐ⱼⱼ, 𝐐ⱼₕ, 𝐐ₕⱼ, 𝐐ₕₕ, 𝐉ᵀ𝐐, 𝐇ᵀ𝐐 = _iitm_project_q_block_values(k, 𝐉, 𝐇, 𝐐)

    ∂𝐐ⱼⱼ, ∂𝐐ⱼₕ,
    ∂𝐐ₕⱼ,
    ∂𝐐ₕₕ = _iitm_project_q_block_derivatives(k, ∂k, 𝐉, 𝐇, 𝐐, ∂𝐉, ∂𝐇, ∂𝐐,
        𝐐ⱼⱼ, 𝐐ⱼₕ, 𝐐ₕⱼ, 𝐐ₕₕ,
        𝐉ᵀ𝐐, 𝐇ᵀ𝐐)

    return 𝐐ⱼⱼ, 𝐐ⱼₕ, 𝐐ₕⱼ, 𝐐ₕₕ, ∂𝐐ⱼⱼ, ∂𝐐ⱼₕ, ∂𝐐ₕⱼ, ∂𝐐ₕₕ
end

function _iitm_project_q_block_derivatives(k, ∂k, 𝐉, 𝐇, 𝐐, ∂𝐉, ∂𝐇, ∂𝐐,
        𝐐ⱼⱼ, 𝐐ⱼₕ, 𝐐ₕⱼ, 𝐐ₕₕ,
        𝐉ᵀ𝐐, 𝐇ᵀ𝐐)
    ∂𝐉ᵀ𝐐 = transpose(∂𝐉) * 𝐐
    ∂𝐇ᵀ𝐐 = transpose(∂𝐇) * 𝐐
    𝐉ᵀ∂𝐐 = transpose(𝐉) * ∂𝐐
    𝐇ᵀ∂𝐐 = transpose(𝐇) * ∂𝐐
    ∂k_over_k = ∂k / k

    ∂𝐐ⱼⱼ = ∂k_over_k * 𝐐ⱼⱼ +
           im * k * (∂𝐉ᵀ𝐐 * 𝐉 + 𝐉ᵀ∂𝐐 * 𝐉 + 𝐉ᵀ𝐐 * ∂𝐉)
    ∂𝐐ⱼₕ = ∂k_over_k * 𝐐ⱼₕ +
           im * k * (∂𝐉ᵀ𝐐 * 𝐇 + 𝐉ᵀ∂𝐐 * 𝐇 + 𝐉ᵀ𝐐 * ∂𝐇)
    ∂𝐐ₕⱼ = ∂k_over_k * 𝐐ₕⱼ +
           im * k * (∂𝐇ᵀ𝐐 * 𝐉 + 𝐇ᵀ∂𝐐 * 𝐉 + 𝐇ᵀ𝐐 * ∂𝐉)
    ∂𝐐ₕₕ = ∂k_over_k * 𝐐ₕₕ +
           im * k * (∂𝐇ᵀ𝐐 * 𝐇 + 𝐇ᵀ∂𝐐 * 𝐇 + 𝐇ᵀ𝐐 * ∂𝐇)

    return ∂𝐐ⱼⱼ, ∂𝐐ⱼₕ, ∂𝐐ₕⱼ, ∂𝐐ₕₕ
end

function _iitm_project_q_block_material_derivatives(k, 𝐉, 𝐇, ∂𝐐)
    𝐉ᵀ∂𝐐 = transpose(𝐉) * ∂𝐐
    𝐇ᵀ∂𝐐 = transpose(𝐇) * ∂𝐐
    ∂𝐐ⱼⱼ = im * k * 𝐉ᵀ∂𝐐 * 𝐉
    ∂𝐐ⱼₕ = im * k * 𝐉ᵀ∂𝐐 * 𝐇
    ∂𝐐ₕⱼ = im * k * 𝐇ᵀ∂𝐐 * 𝐉
    ∂𝐐ₕₕ = im * k * 𝐇ᵀ∂𝐐 * 𝐇

    return ∂𝐐ⱼⱼ, ∂𝐐ⱼₕ, ∂𝐐ₕⱼ, ∂𝐐ₕₕ
end

function _iitm_update_transition_block(𝐓, ∂𝐓, 𝐐ⱼⱼ, 𝐐ⱼₕ, 𝐐ₕⱼ, 𝐐ₕₕ,
        ∂𝐐ⱼⱼ, ∂𝐐ⱼₕ, ∂𝐐ₕⱼ, ∂𝐐ₕₕ)
    return _iitm_update_transition_solve_block(𝐓, ∂𝐓, 𝐐ⱼⱼ, 𝐐ⱼₕ, 𝐐ₕⱼ,
        𝐐ₕₕ, ∂𝐐ⱼⱼ, ∂𝐐ⱼₕ, ∂𝐐ₕⱼ,
        ∂𝐐ₕₕ)
end

function _iitm_transition_solve_cache(𝐓, 𝐐ⱼⱼ, 𝐐ⱼₕ, 𝐐ₕⱼ, 𝐐ₕₕ)
    𝐀 = 𝐈 + 𝐐ⱼₕ
    𝐁 = 𝐈 - 𝐓 * 𝐐ₕₕ
    𝐂 = 𝐈 + 𝐐ₕⱼ
    𝐁_factor = lu(𝐁)
    𝐗 = 𝐁_factor \ (𝐓 * 𝐂)
    𝐓next = 𝐐ⱼⱼ + 𝐀 * 𝐗

    return 𝐓next, 𝐀, 𝐂, 𝐁_factor, 𝐗
end

function _iitm_update_transition_solve_derivative_block(𝐓, ∂𝐓, 𝐐ₕₕ,
        ∂𝐐ⱼⱼ, ∂𝐐ⱼₕ,
        ∂𝐐ₕⱼ, ∂𝐐ₕₕ,
        𝐀, 𝐂, 𝐁_factor,
        𝐗)
    ∂𝐁 = -∂𝐓 * 𝐐ₕₕ - 𝐓 * ∂𝐐ₕₕ
    ∂𝐗 = 𝐁_factor \ (∂𝐓 * 𝐂 + 𝐓 * ∂𝐐ₕⱼ - ∂𝐁 * 𝐗)
    return ∂𝐐ⱼⱼ + ∂𝐐ⱼₕ * 𝐗 + 𝐀 * ∂𝐗
end

function _iitm_update_transition_solve_block(𝐓, ∂𝐓, 𝐐ⱼⱼ, 𝐐ⱼₕ, 𝐐ₕⱼ,
        𝐐ₕₕ, ∂𝐐ⱼⱼ, ∂𝐐ⱼₕ, ∂𝐐ₕⱼ,
        ∂𝐐ₕₕ)
    𝐓next, 𝐀, 𝐂, 𝐁_factor, 𝐗 = _iitm_transition_solve_cache(𝐓, 𝐐ⱼⱼ, 𝐐ⱼₕ, 𝐐ₕⱼ, 𝐐ₕₕ)
    ∂𝐓next = _iitm_update_transition_solve_derivative_block(𝐓, ∂𝐓,
        𝐐ₕₕ, ∂𝐐ⱼⱼ,
        ∂𝐐ⱼₕ, ∂𝐐ₕⱼ,
        ∂𝐐ₕₕ, 𝐀, 𝐂,
        𝐁_factor, 𝐗)
    return 𝐓next, ∂𝐓next
end

function _iitm_repack_axisymmetric_blocks(Ts, nₘₐₓ)
    packed = similar(Ts)
    for m in 0:nₘₐₓ
        nₘᵢₙ = max(1, m)
        nn = nₘₐₓ - nₘᵢₙ + 1
        T′ = similar(Ts[m + 1])

        for n in 1:nn
            for n′ in 1:nn
                T′[n, n′] = Ts[m + 1][2n - 1, 2n′ - 1]
                T′[n, nn + n′] = Ts[m + 1][2n - 1, 2n′]
                T′[nn + n, n′] = Ts[m + 1][2n, 2n′ - 1]
                T′[nn + n, nn + n′] = Ts[m + 1][2n, 2n′]
            end
        end
        packed[m + 1] = T′
    end

    return packed
end

function _iitm_axisymmetric_fixed_geometry_linearization(input, variables)
    s = input.shape
    nₘₐₓ = input.nₘₐₓ
    Nr = input.Nr
    Nϑ_config = input.Nϑ
    p = length(variables)
    RT, CT, k, m, ∂k, ∂m = _iitm_fixed_geometry_parameter_derivatives(input,
        variables)
    wavenumber_indices = _iitm_nonzero_derivative_indices(∂k)
    material_indices = _iitm_nonzero_derivative_indices(∂m)

    rₘᵢₙ = RT(rmin(s))
    rₘₐₓ = RT(rmax(s))
    xr, wr = gausslegendre(RT, Nr)
    @. xr = (rₘₐₓ - rₘᵢₙ) * (xr + 1) / 2 + rₘᵢₙ
    @. wr = (rₘₐₓ - rₘᵢₙ) / 2 * wr

    x, w = gausslegendre(RT, Nϑ_config)
    ϑ = acos.(x)
    Nϑ = has_symmetric_plane(s) ? Nϑ_config ÷ 2 : Nϑ_config
    xᵥ = view(x, 1:Nϑ)
    wᵥ = view(w, 1:Nϑ)
    ϑᵥ = view(ϑ, 1:Nϑ)

    x₀ = k * rₘᵢₙ
    ∂x₀ = [∂k[j] * rₘᵢₙ for j in 1:p]
    a, b, ∂a, ∂b = _bhmie_linearized(RT, x₀, m, nₘₐₓ, ∂x₀, ∂m)
    Ts = Matrix{CT}[]
    ∂Ts_by_var = [[Matrix{CT}(undef, 0, 0) for _ in 0:nₘₐₓ] for _ in 1:p]
    for m_order in 0:nₘₐₓ
        nn = nₘₐₓ + 1 - max(m_order, 1)
        offset = max(m_order, 1) - 1
        𝐓ₘ = zeros(CT, 2nn, 2nn)
        ∂𝐓ₘs = [zeros(CT, 2nn, 2nn) for _ in 1:p]

        for n in 1:nn
            𝐓ₘ[2n - 1, 2n - 1] = -b[n + offset]
            𝐓ₘ[2n, 2n] = -a[n + offset]
            for j in 1:p
                ∂𝐓ₘs[j][2n - 1, 2n - 1] = -∂b[n + offset, j]
                ∂𝐓ₘs[j][2n, 2n] = -∂a[n + offset, j]
            end
        end

        push!(Ts, 𝐓ₘ)
        for j in 1:p
            ∂Ts_by_var[j][m_order + 1] = ∂𝐓ₘs[j]
        end
    end

    d = OffsetArray(zeros(RT, Nϑ, nₘₐₓ + 1, nₘₐₓ + 1), 1:Nϑ, 0:nₘₐₓ,
        0:nₘₐₓ)
    𝜋 = similar(d)
    τ = similar(d)

    Threads.@threads for (i, m_order) in collect(Iterators.product(1:Nϑ, 0:nₘₐₓ))
        wigner_d_recursion!(view(d, i, m_order:nₘₐₓ, m_order), 0, m_order,
            nₘₐₓ, ϑᵥ[i];
            deriv = view(τ, i, m_order:nₘₐₓ, m_order))

        for n in max(m_order, 1):nₘₐₓ
            𝜋[i, n, m_order] = pi_func(RT, m_order, n, ϑᵥ[i];
                d = d[i, n, m_order])
        end
    end

    a½ = [√(RT(n * (n + 1))) for n in 1:nₘₐₓ]
    A = [√(RT(2n + 1) / (2n * (n + 1))) for n in 1:nₘₐₓ]
    sinϑᵥ = sin.(ϑᵥ)

    𝐉 = zeros(CT, 3nₘₐₓ, 2nₘₐₓ)
    𝐇 = zeros(CT, 3nₘₐₓ, 2nₘₐₓ)
    𝐆 = zeros(CT, 3nₘₐₓ, 3nₘₐₓ)
    ∂𝐉s = [similar(𝐉) for _ in 1:p]
    ∂𝐇s = [similar(𝐇) for _ in 1:p]
    ∂𝐆s = [similar(𝐆) for _ in 1:p]

    nₑₓₜᵣₐ = estimate_ricattibesselj_extra_terms(nₘₐₓ, k * rₘₐₓ)
    ψ = zeros(RT, nₘₐₓ)
    z = zeros(RT, nₘₐₓ + nₑₓₜᵣₐ)
    ψ′ = similar(ψ)
    χ = similar(ψ)
    χ′ = similar(ψ)
    ∂ψs = [zeros(CT, nₘₐₓ) for _ in 1:p]
    ∂ψ′s = [zeros(CT, nₘₐₓ) for _ in 1:p]
    ∂χs = [zeros(CT, nₘₐₓ) for _ in 1:p]
    ∂χ′s = [zeros(CT, nₘₐₓ) for _ in 1:p]
    ∂kr = similar(∂k)
    𝐔_by_m = Matrix{CT}[]
    ∂𝐔s_by_m = Vector{Vector{Matrix{CT}}}()
    ∂Us_by_m = Vector{Vector{SMatrix{3, 3, CT, 9}}}()
    for m_order in 0:nₘₐₓ
        nn = nₘₐₓ - max(1, m_order) + 1
        𝐔 = zeros(CT, 3nn, 3nn)
        push!(𝐔_by_m, 𝐔)
        push!(∂𝐔s_by_m, [similar(𝐔) for _ in 1:p])
        push!(∂Us_by_m, [zero(SMatrix{3, 3, CT}) for _ in 1:p])
    end
    within = Vector{Bool}(undef, length(xᵥ))
    ε = Vector{CT}(undef, length(xᵥ))
    m² = m^2

    for (r, wri) in zip(xr, wr)
        kr = k * r
        nₑₓₜᵣₐ = estimate_ricattibesselj_extra_terms(nₘₐₓ, kr)
        ricattibesselj!(ψ, ψ′, z, nₘₐₓ, nₑₓₜᵣₐ, kr)
        ricattibessely!(χ, χ′, nₘₐₓ, kr)

        for j in wavenumber_indices
            ∂kr[j] = ∂k[j] * r
            _iitm_ricatti_argument_derivatives!(∂ψs[j], ∂ψ′s[j], ψ, ψ′, kr,
                ∂kr[j])
            _iitm_ricatti_argument_derivatives!(∂χs[j], ∂χ′s[j], χ, χ′, kr,
                ∂kr[j])
        end

        _iitm_fill_radial_blocks!(𝐉, 𝐇, 𝐆, ∂𝐉s, ∂𝐇s, ∂𝐆s, ψ, ψ′, χ, χ′,
            ∂ψs, ∂ψ′s, ∂χs, ∂χ′s, kr, ∂kr, k, ∂k,
            a½, nₘₐₓ, wavenumber_indices)

        for i in eachindex(xᵥ)
            inside = (r * sinϑᵥ[i], 0, r * xᵥ[i]) ∈ s
            within[i] = inside
            ε[i] = inside ? m² : one(CT)
        end

        for m_order in 0:nₘₐₓ
            nₘᵢₙ = max(1, m_order)
            𝐔 = 𝐔_by_m[m_order + 1]
            ∂𝐔s = ∂𝐔s_by_m[m_order + 1]
            _iitm_axisymmetric_shell_u!(𝐔, ∂𝐔s, s, r, xᵥ, wᵥ, ε, within, d,
                𝜋, τ, a½, A, k, ∂k, m, ∂m, nₘₐₓ,
                m_order, ∂Us_by_m[m_order + 1],
                material_indices)
            𝐉ᵥ = view(𝐉, (3nₘᵢₙ - 2):(3nₘₐₓ), (2nₘᵢₙ - 1):(2nₘₐₓ))
            𝐇ᵥ = view(𝐇, (3nₘᵢₙ - 2):(3nₘₐₓ), (2nₘᵢₙ - 1):(2nₘₐₓ))
            𝐆ᵥ = view(𝐆, (3nₘᵢₙ - 2):(3nₘₐₓ), (3nₘᵢₙ - 2):(3nₘₐₓ))
            𝐑 = 𝐈 - wri * 𝐔 * 𝐆ᵥ
            𝐑_factor = lu(𝐑)
            𝐐 = wri * (𝐑_factor \ 𝐔)

            𝐐ⱼⱼ, 𝐐ⱼₕ, 𝐐ₕⱼ, 𝐐ₕₕ, 𝐉ᵀ𝐐, 𝐇ᵀ𝐐 = _iitm_project_q_block_values(k, 𝐉ᵥ, 𝐇ᵥ, 𝐐)

            𝐓_old = Ts[m_order + 1]
            𝐓_next, 𝐀, 𝐂, 𝐁_factor,
            𝐗 = _iitm_transition_solve_cache(𝐓_old, 𝐐ⱼⱼ, 𝐐ⱼₕ, 𝐐ₕⱼ,
                𝐐ₕₕ)
            Ts[m_order + 1] = 𝐓_next

            for j in 1:p
                if iszero(∂k[j])
                    ∂𝐑 = -wri * (∂𝐔s[j] * 𝐆ᵥ)
                    ∂𝐐 = 𝐑_factor \ (wri * ∂𝐔s[j] - ∂𝐑 * 𝐐)
                    ∂𝐐ⱼⱼ, ∂𝐐ⱼₕ,
                    ∂𝐐ₕⱼ, ∂𝐐ₕₕ = _iitm_project_q_block_material_derivatives(k, 𝐉ᵥ,
                        𝐇ᵥ, ∂𝐐)
                else
                    ∂𝐉ᵥ = view(∂𝐉s[j], (3nₘᵢₙ - 2):(3nₘₐₓ),
                        (2nₘᵢₙ - 1):(2nₘₐₓ))
                    ∂𝐇ᵥ = view(∂𝐇s[j], (3nₘᵢₙ - 2):(3nₘₐₓ),
                        (2nₘᵢₙ - 1):(2nₘₐₓ))
                    ∂𝐆ᵥ = view(∂𝐆s[j], (3nₘᵢₙ - 2):(3nₘₐₓ),
                        (3nₘᵢₙ - 2):(3nₘₐₓ))
                    ∂𝐑 = -wri * (∂𝐔s[j] * 𝐆ᵥ + 𝐔 * ∂𝐆ᵥ)
                    ∂𝐐 = 𝐑_factor \ (wri * ∂𝐔s[j] - ∂𝐑 * 𝐐)
                    ∂𝐐ⱼⱼ, ∂𝐐ⱼₕ,
                    ∂𝐐ₕⱼ,
                    ∂𝐐ₕₕ = _iitm_project_q_block_derivatives(k, ∂k[j], 𝐉ᵥ, 𝐇ᵥ,
                        𝐐, ∂𝐉ᵥ, ∂𝐇ᵥ,
                        ∂𝐐, 𝐐ⱼⱼ, 𝐐ⱼₕ,
                        𝐐ₕⱼ, 𝐐ₕₕ, 𝐉ᵀ𝐐,
                        𝐇ᵀ𝐐)
                end
                ∂𝐓_next = _iitm_update_transition_solve_derivative_block(𝐓_old,
                    ∂Ts_by_var[j][m_order + 1],
                    𝐐ₕₕ,
                    ∂𝐐ⱼⱼ,
                    ∂𝐐ⱼₕ,
                    ∂𝐐ₕⱼ,
                    ∂𝐐ₕₕ,
                    𝐀, 𝐂,
                    𝐁_factor,
                    𝐗)
                ∂Ts_by_var[j][m_order + 1] = ∂𝐓_next
            end
        end
    end

    value_blocks = _iitm_repack_axisymmetric_blocks(Ts, nₘₐₓ)
    jacobian = map(∂Ts_by_var) do ∂Ts
        derivative_blocks = _iitm_repack_axisymmetric_blocks(∂Ts, nₘₐₓ)
        AxisymmetricTransitionMatrix{CT, nₘₐₓ, typeof(derivative_blocks), RT}(derivative_blocks)
    end
    value = AxisymmetricTransitionMatrix{CT, nₘₐₓ, typeof(value_blocks), RT}(value_blocks)

    return LinearizationResult(value, jacobian, variables;
        metadata = (; backend = :iitm_axisymmetric_analytic,
            variant = :axisymmetric,
            λ = input.λ,
            nₘₐₓ,
            Nr,
            Nϑ = Nϑ_config))
end

_iitm_nfold_period(::AbstractNFoldShape{N}) where {N} = N

function _iitm_effective_variant(input, backend::IITMLinearization)
    backend.variant != :auto && return backend.variant
    input.shape isa AbstractAxisymmetricShape && return :axisymmetric
    input.shape isa AbstractNFoldShape && return :nfold
    return :arbitrary
end

function _iitm_fill_ordered_radial_blocks!(𝐉, 𝐇, 𝐆, ∂𝐉s, ∂𝐇s, ∂𝐆s,
        ψ, ψ′, χ, χ′, ∂ψs, ∂ψ′s,
        ∂χs, ∂χ′s, kr, ∂kr, k, ∂k,
        a½, order_degree,
        wavenumber_indices = eachindex(∂k))
    CT = eltype(𝐉)
    fill!(𝐉, zero(CT))
    fill!(𝐇, zero(CT))
    fill!(𝐆, zero(CT))
    for ∂𝐉 in ∂𝐉s
        fill!(∂𝐉, zero(CT))
    end
    for ∂𝐇 in ∂𝐇s
        fill!(∂𝐇, zero(CT))
    end
    for ∂𝐆 in ∂𝐆s
        fill!(∂𝐆, zero(CT))
    end

    for (i, (n, _)) in order_degree
        𝐉ᵈ = @SMatrix [ψ[n]/kr 0
                       0 ψ′[n]/kr
                       0 a½[n] * ψ[n]/kr^2]
        𝐘ᵈ = @SMatrix [χ[n]/kr 0
                       0 χ′[n]/kr
                       0 a½[n] * χ[n]/kr^2]
        𝐇ᵈ = 𝐉ᵈ + 1im * 𝐘ᵈ
        𝐆ᵈ_base = 𝐇ᵈ * transpose(𝐉ᵈ) + 𝐉ᵈ * transpose(𝐇ᵈ)

        𝐉[(3i - 2):(3i), (2i - 1):(2i)] .= 𝐉ᵈ
        𝐇[(3i - 2):(3i), (2i - 1):(2i)] .= 𝐇ᵈ
        𝐆[(3i - 2):(3i), (3i - 2):(3i)] .= 𝐆ᵈ_base .* (im * k / 2)

        for j in wavenumber_indices
            ∂krⱼ = ∂kr[j]
            ∂𝐉ᵈ = @SMatrix [∂ψs[j][n]/kr - ψ[n] * ∂krⱼ/kr^2 0
                            0 ∂ψ′s[j][n]/kr - ψ′[n] * ∂krⱼ/kr^2
                            0 a½[n] * (∂ψs[j][n]/kr^2 - 2ψ[n] * ∂krⱼ/kr^3)]
            ∂𝐘ᵈ = @SMatrix [∂χs[j][n]/kr - χ[n] * ∂krⱼ/kr^2 0
                            0 ∂χ′s[j][n]/kr - χ′[n] * ∂krⱼ/kr^2
                            0 a½[n] * (∂χs[j][n]/kr^2 - 2χ[n] * ∂krⱼ/kr^3)]
            ∂𝐇ᵈ = ∂𝐉ᵈ + 1im * ∂𝐘ᵈ
            ∂𝐆ᵈ_base = ∂𝐇ᵈ * transpose(𝐉ᵈ) + 𝐇ᵈ * transpose(∂𝐉ᵈ) +
                       ∂𝐉ᵈ * transpose(𝐇ᵈ) + 𝐉ᵈ * transpose(∂𝐇ᵈ)
            ∂𝐆ᵈ = ∂𝐆ᵈ_base .* (im * k / 2) .+
                  𝐆ᵈ_base .* (im * ∂k[j] / 2)

            ∂𝐉s[j][(3i - 2):(3i), (2i - 1):(2i)] .= ∂𝐉ᵈ
            ∂𝐇s[j][(3i - 2):(3i), (2i - 1):(2i)] .= ∂𝐇ᵈ
            ∂𝐆s[j][(3i - 2):(3i), (3i - 2):(3i)] .= ∂𝐆ᵈ
        end
    end

    return nothing
end

function _iitm_material_tables(s, r, x, ϑ, xφ, m, ∂m, material_indices)
    CT = typeof(m)
    ε = Matrix{CT}(undef, length(xφ), length(x))
    ∂εs = [Matrix{CT}(undef, length(xφ), length(x)) for _ in material_indices]

    for (jφ, φ) in enumerate(xφ), i in eachindex(x)

        local_m = refractive_index(s,
            (r * sin(ϑ[i]) * cos(φ),
                r * sin(ϑ[i]) * sin(φ),
                r * x[i]))
        ε[jφ, i] = local_m^2
        material_matches = local_m == m
        for (jm, j) in enumerate(material_indices)
            ∂local_m = material_matches ? ∂m[j] : zero(CT)
            ∂εs[jm][jφ, i] = 2local_m * ∂local_m
        end
    end

    return ε, ∂εs
end

function _iitm_collect_fourier_coefficients(spectrum, nₘₐₓ, wφ, mode_bins)
    _, Nϑ = size(spectrum)
    qs, bins = mode_bins
    coeff = OffsetArray(zeros(ComplexF64, 4nₘₐₓ + 1, Nϑ),
        (-2nₘₐₓ):(2nₘₐₓ), 1:Nϑ)

    for i in 1:Nϑ
        for iq in eachindex(qs)
            q = qs[iq]
            idx = mod(-bins[iq], size(spectrum, 1)) + 1
            coeff[q, i] = wφ * spectrum[idx, i]
        end
    end

    return coeff
end

function _iitm_fourier_coefficients_with_derivatives(ε, ∂εs, nₘₐₓ, wφ,
        workspace,
        mode_bins)
    coeff_ε, coeff_εinv = _azimuthal_fourier_coefficients(ε, nₘₐₓ, wφ,
        workspace,
        mode_bins)
    ∂coeffs = map(∂εs) do ∂ε
        @. workspace.contrast = ∂ε
        @. workspace.contrast_inv = ∂ε / ε^2
        mul!(workspace.spectrum, workspace.plan, workspace.contrast)
        mul!(workspace.spectrum_inv, workspace.plan, workspace.contrast_inv)
        (_iitm_collect_fourier_coefficients(workspace.spectrum, nₘₐₓ, wφ,
                mode_bins),
            _iitm_collect_fourier_coefficients(workspace.spectrum_inv, nₘₐₓ, wφ,
                mode_bins))
    end

    return coeff_ε, coeff_εinv, ∂coeffs
end

function _iitm_ordered_shell_data(s, r, x, ϑ, xφ, m, ∂m, nₘₐₓ, wφ,
        fourier_workspace, fourier_modes,
        material_indices = _iitm_nonzero_derivative_indices(∂m))
    ε, ∂εs = _iitm_material_tables(s, r, x, ϑ, xφ, m, ∂m,
        material_indices)
    fourier_coeffs = isnothing(fourier_workspace) ? nothing :
                     _iitm_fourier_coefficients_with_derivatives(ε, ∂εs,
        nₘₐₓ, wφ,
        fourier_workspace,
        fourier_modes)
    return ε, ∂εs, fourier_coeffs, material_indices
end

function _iitm_ordered_shell_u!(𝐔, ∂𝐔s, s, r, x, w, ϑ, xφ, wφ, d, 𝜋, τ,
        a½, A, k, ∂k, m, ∂m, nₘₐₓ, order_degree,
        scale_factor, fourier_workspace,
        fourier_modes, shell_data = nothing,
        ∂Us = nothing)
    CT = typeof(m)
    pvars = length(∂m)
    fill!(𝐔, zero(CT))
    for ∂𝐔 in ∂𝐔s
        fill!(∂𝐔, zero(CT))
    end
    kr = k * r
    ε, ∂εs,
    fourier_coeffs,
    material_indices = isnothing(shell_data) ?
                       _iitm_ordered_shell_data(s, r, x, ϑ, xφ, m, ∂m, nₘₐₓ, wφ,
        fourier_workspace, fourier_modes) :
                       shell_data
    if isnothing(∂Us)
        ∂Us = [zero(SMatrix{3, 3, CT}) for _ in 1:pvars]
    end

    for (q, (n′, m′)) in order_degree
        for (p, (n, m_order)) in order_degree
            sig = iseven(m_order + m′) ? 1 : -1
            if has_symmetric_plane(s)
                c = iseven(n + m_order + n′ + m′) ? 2 : 0
                c̃ = 2 - c
            else
                c = 1
                c̃ = 1
            end

            U = zero(SMatrix{3, 3, CT})
            fill!(∂Us, zero(SMatrix{3, 3, CT}))
            freq = m′ - m_order

            for i in eachindex(x)
                pptt = 𝜋[i, n, m_order] * 𝜋[i, n′, m′] +
                       τ[i, n, m_order] * τ[i, n′, m′]
                pttp = 𝜋[i, n, m_order] * τ[i, n′, m′] +
                       τ[i, n, m_order] * 𝜋[i, n′, m′]
                dd = d[i, n, m_order] * d[i, n′, m′]

                if isnothing(fourier_coeffs)
                    for (jφ, φ) in enumerate(xφ)
                        wφⱼ = wφ isa Number ? wφ : wφ[jφ]
                        phase = cis(freq * φ)
                        ΔU = @SMatrix [c*pptt -c̃*im*pttp 0
                                       c̃*im*pttp c*pptt 0
                                       0 0 c * a½[n] * a½[n′] * dd/ε[jφ, i]]
                        contrast = ε[jφ, i] - 1
                        U += w[i] * wφⱼ * phase * contrast * ΔU
                        for (jm, j) in enumerate(material_indices)
                            ∂ε = ∂εs[jm][jφ, i]
                            ∂ΔU = @SMatrix [zero(CT) zero(CT) zero(CT)
                                            zero(CT) zero(CT) zero(CT)
                                            zero(CT) zero(CT) -c * a½[n] * a½[n′] * dd * ∂ε / ε[jφ, i]^2]
                            ∂Us[j] += w[i] * wφⱼ * phase *
                                      (∂ε * ΔU + contrast * ∂ΔU)
                        end
                    end
                else
                    coeff_ε, coeff_εinv, ∂coeffs = fourier_coeffs
                    cε = coeff_ε[freq, i]
                    cεinv = coeff_εinv[freq, i]
                    U += w[i] * @SMatrix [c*pptt*cε -c̃*im*pttp*cε 0
                                   c̃*im*pttp*cε c*pptt*cε 0
                                   0 0 c * a½[n] * a½[n′] * dd * cεinv]
                    for (jm, j) in enumerate(material_indices)
                        ∂cε, ∂cεinv = ∂coeffs[jm]
                        dcε = ∂cε[freq, i]
                        dcεinv = ∂cεinv[freq, i]
                        ∂Us[j] += w[i] * @SMatrix [c*pptt*dcε -c̃*im*pttp*dcε 0
                                            c̃*im*pttp*dcε c*pptt*dcε 0
                                            0 0 c * a½[n] * a½[n′] * dd * dcεinv]
                    end
                end
            end

            prefactor = kr^2 * sig * A[n] * A[n′] * scale_factor
            U_unscaled = U
            U = prefactor * U_unscaled
            for j in 1:pvars
                ∂prefactor = 2kr * ∂k[j] * r * sig * A[n] * A[n′] *
                             scale_factor
                ∂Us[j] = ∂prefactor * U_unscaled + prefactor * ∂Us[j]
            end

            rows = (3p - 2):(3p)
            cols = (3q - 2):(3q)
            𝐔[rows, cols] .= U
            for j in 1:pvars
                ∂𝐔s[j][rows, cols] .= ∂Us[j]
            end
        end
    end

    return 𝐔, ∂𝐔s
end

function _iitm_ordered_shell_u(s, r, x, w, ϑ, xφ, wφ, d, 𝜋, τ, a½, A,
        k, ∂k, m, ∂m, nₘₐₓ, order_degree,
        scale_factor, fourier_workspace,
        fourier_modes)
    CT = typeof(m)
    L = length(order_degree)
    𝐔 = zeros(CT, 3L, 3L)
    ∂𝐔s = [zeros(CT, 3L, 3L) for _ in eachindex(∂m)]
    return _iitm_ordered_shell_u!(𝐔, ∂𝐔s, s, r, x, w, ϑ, xφ, wφ, d, 𝜋, τ,
        a½, A, k, ∂k, m, ∂m, nₘₐₓ, order_degree,
        scale_factor, fourier_workspace,
        fourier_modes)
end

function _iitm_update_general_block!(𝐓, ∂𝐓s, 𝐉, 𝐇, 𝐆, 𝐔, ∂𝐉s, ∂𝐇s,
        ∂𝐆s, ∂𝐔s, wri, k, ∂k,
        𝐓_old = copy(𝐓))
    𝐑 = 𝐈 - wri * 𝐔 * 𝐆
    𝐑_factor = lu(𝐑)
    𝐐 = wri * (𝐑_factor \ 𝐔)
    𝐐ⱼⱼ, 𝐐ⱼₕ, 𝐐ₕⱼ, 𝐐ₕₕ, 𝐉ᵀ𝐐, 𝐇ᵀ𝐐 = _iitm_project_q_block_values(k, 𝐉, 𝐇, 𝐐)
    copyto!(𝐓_old, 𝐓)
    𝐓_next, 𝐀, 𝐂, 𝐁_factor, 𝐗 = _iitm_transition_solve_cache(𝐓_old, 𝐐ⱼⱼ, 𝐐ⱼₕ, 𝐐ₕⱼ, 𝐐ₕₕ)
    𝐓 .= 𝐓_next

    for j in eachindex(∂𝐓s)
        ∂𝐑 = iszero(∂k[j]) ? -wri * (∂𝐔s[j] * 𝐆) :
             -wri * (∂𝐔s[j] * 𝐆 + 𝐔 * ∂𝐆s[j])
        ∂𝐐 = 𝐑_factor \ (wri * ∂𝐔s[j] - ∂𝐑 * 𝐐)
        if iszero(∂k[j])
            ∂𝐐ⱼⱼ, ∂𝐐ⱼₕ, ∂𝐐ₕⱼ, ∂𝐐ₕₕ = _iitm_project_q_block_material_derivatives(k, 𝐉, 𝐇, ∂𝐐)
        else
            ∂𝐐ⱼⱼ, ∂𝐐ⱼₕ,
            ∂𝐐ₕⱼ,
            ∂𝐐ₕₕ = _iitm_project_q_block_derivatives(k, ∂k[j], 𝐉, 𝐇, 𝐐, ∂𝐉s[j],
                ∂𝐇s[j], ∂𝐐, 𝐐ⱼⱼ, 𝐐ⱼₕ,
                𝐐ₕⱼ, 𝐐ₕₕ, 𝐉ᵀ𝐐, 𝐇ᵀ𝐐)
        end
        ∂𝐓s[j] .= _iitm_update_transition_solve_derivative_block(𝐓_old, ∂𝐓s[j],
            𝐐ₕₕ, ∂𝐐ⱼⱼ,
            ∂𝐐ⱼₕ, ∂𝐐ₕⱼ,
            ∂𝐐ₕₕ, 𝐀, 𝐂,
            𝐁_factor, 𝐗)
    end

    return nothing
end

function _iitm_repack_transition_matrix(𝐓, order_degree, nₘₐₓ)
    CT = eltype(𝐓)
    𝐓′ = OffsetArray(zeros(CT, 2nₘₐₓ + 1, nₘₐₓ, 2nₘₐₓ + 1, nₘₐₓ, 2, 2),
        (-nₘₐₓ):nₘₐₓ, 1:nₘₐₓ, (-nₘₐₓ):nₘₐₓ, 1:nₘₐₓ, 1:2,
        1:2)

    for (j, (n′, m′)) in order_degree
        for (i, (n, m_order)) in order_degree
            𝐓′[m_order, n, m′, n′, 1, 1] = 𝐓[2i - 1, 2j - 1]
            𝐓′[m_order, n, m′, n′, 1, 2] = 𝐓[2i - 1, 2j]
            𝐓′[m_order, n, m′, n′, 2, 1] = 𝐓[2i, 2j - 1]
            𝐓′[m_order, n, m′, n′, 2, 2] = 𝐓[2i, 2j]
        end
    end

    return 𝐓′
end

function _iitm_arbitrary_fixed_geometry_linearization(input, variables)
    s = input.shape
    nₘₐₓ = input.nₘₐₓ
    Nr = input.Nr
    Nϑ_config = input.Nϑ
    Nφ = input.Nφ
    pvars = length(variables)
    RT, CT, k, m, ∂k, ∂m = _iitm_fixed_geometry_parameter_derivatives(input,
        variables)
    wavenumber_indices = _iitm_nonzero_derivative_indices(∂k)
    material_indices = _iitm_nonzero_derivative_indices(∂m)

    rₘᵢₙ = RT(rmin(s))
    rₘₐₓ = RT(rmax(s))
    xr, wr = gausslegendre(RT, Nr)
    @. xr = (rₘₐₓ - rₘᵢₙ) * (xr + 1) / 2 + rₘᵢₙ
    @. wr = (rₘₐₓ - rₘᵢₙ) / 2 * wr

    x, w = gausslegendre(RT, Nϑ_config)
    ϑ = acos.(x)
    Nϑ = has_symmetric_plane(s) ? Nϑ_config ÷ 2 : Nϑ_config
    xᵥ = view(x, 1:Nϑ)
    wᵥ = view(w, 1:Nϑ)
    ϑᵥ = view(ϑ, 1:Nϑ)

    xφ = range(zero(RT), 2 * RT(π), length = Nφ + 1)[1:(end - 1)]
    wφ = 2 * RT(π) / Nφ
    fourier_workspace = CT <: ComplexF64 ? _azimuthal_fourier_workspace(Nφ, Nϑ, nₘₐₓ) :
                        nothing
    fourier_modes = CT <: ComplexF64 ? _azimuthal_fourier_mode_bins(nₘₐₓ) :
                    nothing

    order_degree = collect(enumerate(collect(OrderDegreeIterator(nₘₐₓ))))
    L = length(order_degree)
    𝐓 = zeros(CT, 2L, 2L)
    ∂𝐓s = [zeros(CT, 2L, 2L) for _ in 1:pvars]
    ∂x₀ = [∂k[j] * rₘᵢₙ for j in 1:pvars]
    a, b, ∂a, ∂b = _bhmie_linearized(RT, k * rₘᵢₙ, m, nₘₐₓ, ∂x₀, ∂m)
    for (i, (n, _)) in order_degree
        𝐓[2i - 1, 2i - 1] = -b[n]
        𝐓[2i, 2i] = -a[n]
        for j in 1:pvars
            ∂𝐓s[j][2i - 1, 2i - 1] = -∂b[n, j]
            ∂𝐓s[j][2i, 2i] = -∂a[n, j]
        end
    end

    d = OffsetArray(zeros(RT, Nϑ, nₘₐₓ + 1, 2nₘₐₓ + 1), 1:Nϑ, 0:nₘₐₓ,
        (-nₘₐₓ):nₘₐₓ)
    𝜋 = similar(d)
    τ = similar(d)
    Threads.@threads for (i, m_order) in collect(Iterators.product(1:Nϑ,
        (-nₘₐₓ):nₘₐₓ))
        wigner_d_recursion!(view(d, i, abs(m_order):nₘₐₓ, m_order), 0,
            m_order, nₘₐₓ, ϑᵥ[i];
            deriv = view(τ, i, abs(m_order):nₘₐₓ, m_order))
        for n in max(abs(m_order), 1):nₘₐₓ
            𝜋[i, n, m_order] = pi_func(RT, m_order, n, ϑᵥ[i];
                d = d[i, n, m_order])
        end
    end

    a½ = [√(RT(n * (n + 1))) for n in 1:nₘₐₓ]
    A = [√(RT(2n + 1) / (2n * (n + 1))) for n in 1:nₘₐₓ]
    𝐉 = zeros(CT, 3L, 2L)
    𝐇 = zeros(CT, 3L, 2L)
    𝐆 = zeros(CT, 3L, 3L)
    𝐔 = zeros(CT, 3L, 3L)
    ∂𝐉s = [similar(𝐉) for _ in 1:pvars]
    ∂𝐇s = [similar(𝐇) for _ in 1:pvars]
    ∂𝐆s = [similar(𝐆) for _ in 1:pvars]
    ∂𝐔s = [similar(𝐔) for _ in 1:pvars]
    ψ = zeros(RT, nₘₐₓ)
    z = zeros(RT, nₘₐₓ + estimate_ricattibesselj_extra_terms(nₘₐₓ, k * rₘₐₓ))
    ψ′ = similar(ψ)
    χ = similar(ψ)
    χ′ = similar(ψ)
    ∂ψs = [zeros(CT, nₘₐₓ) for _ in 1:pvars]
    ∂ψ′s = [zeros(CT, nₘₐₓ) for _ in 1:pvars]
    ∂χs = [zeros(CT, nₘₐₓ) for _ in 1:pvars]
    ∂χ′s = [zeros(CT, nₘₐₓ) for _ in 1:pvars]
    ∂kr = similar(∂k)
    𝐓_old = similar(𝐓)
    ∂Us = [zero(SMatrix{3, 3, CT}) for _ in 1:pvars]

    for (r, wri) in zip(xr, wr)
        kr = k * r
        nₑₓₜᵣₐ = estimate_ricattibesselj_extra_terms(nₘₐₓ, kr)
        length(z) < nₘₐₓ + nₑₓₜᵣₐ && resize!(z, nₘₐₓ + nₑₓₜᵣₐ)
        ricattibesselj!(ψ, ψ′, z, nₘₐₓ, nₑₓₜᵣₐ, kr)
        ricattibessely!(χ, χ′, nₘₐₓ, kr)
        for j in wavenumber_indices
            ∂kr[j] = ∂k[j] * r
            _iitm_ricatti_argument_derivatives!(∂ψs[j], ∂ψ′s[j], ψ, ψ′, kr,
                ∂kr[j])
            _iitm_ricatti_argument_derivatives!(∂χs[j], ∂χ′s[j], χ, χ′, kr,
                ∂kr[j])
        end
        _iitm_fill_ordered_radial_blocks!(𝐉, 𝐇, 𝐆, ∂𝐉s, ∂𝐇s, ∂𝐆s, ψ, ψ′,
            χ, χ′, ∂ψs, ∂ψ′s, ∂χs, ∂χ′s, kr,
            ∂kr, k, ∂k, a½, order_degree,
            wavenumber_indices)
        shell_data = _iitm_ordered_shell_data(s, r, xᵥ, ϑᵥ, xφ, m, ∂m,
            nₘₐₓ, wφ, fourier_workspace,
            fourier_modes, material_indices)
        _iitm_ordered_shell_u!(𝐔, ∂𝐔s, s, r, xᵥ, wᵥ, ϑᵥ, xφ, wφ, d, 𝜋,
            τ, a½, A, k, ∂k, m, ∂m, nₘₐₓ, order_degree,
            inv(2 * RT(π)), fourier_workspace,
            fourier_modes, shell_data, ∂Us)
        _iitm_update_general_block!(𝐓, ∂𝐓s, 𝐉, 𝐇, 𝐆, 𝐔, ∂𝐉s, ∂𝐇s, ∂𝐆s,
            ∂𝐔s, wri, k, ∂k, 𝐓_old)
    end

    value_container = _iitm_repack_transition_matrix(𝐓, order_degree, nₘₐₓ)
    value = TransitionMatrix{CT, nₘₐₓ, typeof(value_container)}(value_container)
    jacobian = map(∂𝐓s) do ∂𝐓
        container = _iitm_repack_transition_matrix(∂𝐓, order_degree, nₘₐₓ)
        TransitionMatrix{CT, nₘₐₓ, typeof(container)}(container)
    end

    return LinearizationResult(value, jacobian, variables;
        metadata = (; backend = :iitm_arbitrary_analytic,
            variant = :arbitrary,
            λ = input.λ,
            nₘₐₓ,
            Nr,
            Nϑ = Nϑ_config,
            Nφ))
end

function _iitm_nfold_fixed_geometry_linearization(input, variables)
    s = input.shape
    N = _iitm_nfold_period(s)
    nₘₐₓ = input.nₘₐₓ
    Nr = input.Nr
    Nϑ_config = input.Nϑ
    Nφ = input.Nφ
    pvars = length(variables)
    RT, CT, k, m, ∂k, ∂m = _iitm_fixed_geometry_parameter_derivatives(input,
        variables)
    wavenumber_indices = _iitm_nonzero_derivative_indices(∂k)
    material_indices = _iitm_nonzero_derivative_indices(∂m)

    rₘᵢₙ = RT(rmin(s))
    rₘₐₓ = RT(rmax(s))
    xr, wr = gausslegendre(RT, Nr)
    @. xr = (rₘₐₓ - rₘᵢₙ) * (xr + 1) / 2 + rₘᵢₙ
    @. wr = (rₘₐₓ - rₘᵢₙ) / 2 * wr

    x, w = gausslegendre(RT, Nϑ_config)
    ϑ = acos.(x)
    Nϑ = has_symmetric_plane(s) ? Nϑ_config ÷ 2 : Nϑ_config
    xᵥ = view(x, 1:Nϑ)
    wᵥ = view(w, 1:Nϑ)
    ϑᵥ = view(ϑ, 1:Nϑ)

    if CT <: ComplexF64
        xφ = range(zero(RT), 2 * RT(π) / N, length = Nφ + 1)[1:(end - 1)]
        wφ = 2 * RT(π) / (N * Nφ)
    else
        xφ, wφ = gausslegendre(RT, Nφ)
        @. xφ = (xφ + 1) * π / N
        @. wφ = π / N * wφ
    end
    fourier_workspace = CT <: ComplexF64 ? _azimuthal_fourier_workspace(Nφ, Nϑ, nₘₐₓ) :
                        nothing
    fourier_modes = CT <: ComplexF64 ? _azimuthal_fourier_mode_bins(nₘₐₓ, N) :
                    nothing

    all_orders = collect(OrderDegreeIterator(nₘₐₓ))
    order_groups = [collect(enumerate([nm for nm in all_orders
                                       if (nm[2] % N + N) % N == i]))
                    for i in 0:(N - 1)]
    ∂x₀ = [∂k[j] * rₘᵢₙ for j in 1:pvars]
    a, b, ∂a, ∂b = _bhmie_linearized(RT, k * rₘᵢₙ, m, nₘₐₓ, ∂x₀, ∂m)
    𝐓s = [zeros(CT, 2length(group), 2length(group)) for group in order_groups]
    ∂𝐓s_by_var = [[zeros(CT, size(𝐓s[g])...) for g in eachindex(𝐓s)]
                  for _ in 1:pvars]

    for (g, group) in enumerate(order_groups)
        for (i, (n, _)) in group
            𝐓s[g][2i - 1, 2i - 1] = -b[n]
            𝐓s[g][2i, 2i] = -a[n]
            for j in 1:pvars
                ∂𝐓s_by_var[j][g][2i - 1, 2i - 1] = -∂b[n, j]
                ∂𝐓s_by_var[j][g][2i, 2i] = -∂a[n, j]
            end
        end
    end
    ∂𝐓s_by_group = [[∂𝐓s_by_var[j][g] for j in 1:pvars]
                    for g in eachindex(𝐓s)]

    d = OffsetArray(zeros(RT, Nϑ, nₘₐₓ + 1, 2nₘₐₓ + 1), 1:Nϑ, 0:nₘₐₓ,
        (-nₘₐₓ):nₘₐₓ)
    𝜋 = similar(d)
    τ = similar(d)
    Threads.@threads for (i, m_order) in collect(Iterators.product(1:Nϑ,
        (-nₘₐₓ):nₘₐₓ))
        wigner_d_recursion!(view(d, i, abs(m_order):nₘₐₓ, m_order), 0,
            m_order, nₘₐₓ, ϑᵥ[i];
            deriv = view(τ, i, abs(m_order):nₘₐₓ, m_order))
        for n in max(abs(m_order), 1):nₘₐₓ
            𝜋[i, n, m_order] = pi_func(RT, m_order, n, ϑᵥ[i];
                d = d[i, n, m_order])
        end
    end

    a½ = [√(RT(n * (n + 1))) for n in 1:nₘₐₓ]
    A = [√(RT(2n + 1) / (2n * (n + 1))) for n in 1:nₘₐₓ]
    ψ = zeros(RT, nₘₐₓ)
    z = zeros(RT, nₘₐₓ + estimate_ricattibesselj_extra_terms(nₘₐₓ, k * rₘₐₓ))
    ψ′ = similar(ψ)
    χ = similar(ψ)
    χ′ = similar(ψ)
    ∂ψs = [zeros(CT, nₘₐₓ) for _ in 1:pvars]
    ∂ψ′s = [zeros(CT, nₘₐₓ) for _ in 1:pvars]
    ∂χs = [zeros(CT, nₘₐₓ) for _ in 1:pvars]
    ∂χ′s = [zeros(CT, nₘₐₓ) for _ in 1:pvars]
    𝐉_groups = [zeros(CT, 3length(group), 2length(group)) for group in order_groups]
    𝐇_groups = [zeros(CT, 3length(group), 2length(group)) for group in order_groups]
    𝐆_groups = [zeros(CT, 3length(group), 3length(group)) for group in order_groups]
    𝐔_groups = [zeros(CT, 3length(group), 3length(group)) for group in order_groups]
    𝐓_old_groups = [similar(𝐓s[g]) for g in eachindex(𝐓s)]
    ∂𝐉_groups = [[similar(𝐉_groups[g]) for _ in 1:pvars]
                 for g in eachindex(order_groups)]
    ∂𝐇_groups = [[similar(𝐇_groups[g]) for _ in 1:pvars]
                 for g in eachindex(order_groups)]
    ∂𝐆_groups = [[similar(𝐆_groups[g]) for _ in 1:pvars]
                 for g in eachindex(order_groups)]
    ∂𝐔_groups = [[similar(𝐔_groups[g]) for _ in 1:pvars]
                 for g in eachindex(order_groups)]
    ∂Us_groups = [[zero(SMatrix{3, 3, CT}) for _ in 1:pvars]
                  for _ in eachindex(order_groups)]
    ∂kr = similar(∂k)

    for (r, wri) in zip(xr, wr)
        kr = k * r
        nₑₓₜᵣₐ = estimate_ricattibesselj_extra_terms(nₘₐₓ, kr)
        length(z) < nₘₐₓ + nₑₓₜᵣₐ && resize!(z, nₘₐₓ + nₑₓₜᵣₐ)
        ricattibesselj!(ψ, ψ′, z, nₘₐₓ, nₑₓₜᵣₐ, kr)
        ricattibessely!(χ, χ′, nₘₐₓ, kr)
        for j in wavenumber_indices
            ∂kr[j] = ∂k[j] * r
            _iitm_ricatti_argument_derivatives!(∂ψs[j], ∂ψ′s[j], ψ, ψ′, kr,
                ∂kr[j])
            _iitm_ricatti_argument_derivatives!(∂χs[j], ∂χ′s[j], χ, χ′, kr,
                ∂kr[j])
        end
        shell_data = _iitm_ordered_shell_data(s, r, xᵥ, ϑᵥ, xφ, m, ∂m,
            nₘₐₓ, wφ, fourier_workspace,
            fourier_modes, material_indices)

        for (g, group) in enumerate(order_groups)
            𝐉 = 𝐉_groups[g]
            𝐇 = 𝐇_groups[g]
            𝐆 = 𝐆_groups[g]
            𝐔 = 𝐔_groups[g]
            ∂𝐉s = ∂𝐉_groups[g]
            ∂𝐇s = ∂𝐇_groups[g]
            ∂𝐆s = ∂𝐆_groups[g]
            ∂𝐔s = ∂𝐔_groups[g]
            _iitm_fill_ordered_radial_blocks!(𝐉, 𝐇, 𝐆, ∂𝐉s, ∂𝐇s, ∂𝐆s,
                ψ, ψ′, χ, χ′, ∂ψs, ∂ψ′s,
                ∂χs, ∂χ′s, kr, ∂kr, k, ∂k,
                a½, group, wavenumber_indices)
            _iitm_ordered_shell_u!(𝐔, ∂𝐔s, s, r, xᵥ, wᵥ, ϑᵥ, xφ, wφ, d,
                𝜋, τ, a½, A, k, ∂k, m, ∂m, nₘₐₓ, group,
                RT(N) / (2 * RT(π)), fourier_workspace,
                fourier_modes, shell_data, ∂Us_groups[g])
            _iitm_update_general_block!(𝐓s[g], ∂𝐓s_by_group[g],
                𝐉, 𝐇, 𝐆, 𝐔, ∂𝐉s, ∂𝐇s, ∂𝐆s, ∂𝐔s,
                wri, k, ∂k, 𝐓_old_groups[g])
        end
    end

    value_container = OffsetArray(zeros(CT, 2nₘₐₓ + 1, nₘₐₓ, 2nₘₐₓ + 1,
            nₘₐₓ, 2, 2),
        (-nₘₐₓ):nₘₐₓ, 1:nₘₐₓ, (-nₘₐₓ):nₘₐₓ,
        1:nₘₐₓ, 1:2, 1:2)
    derivative_containers = [similar(value_container) for _ in 1:pvars]
    for container in derivative_containers
        fill!(container, zero(CT))
    end

    for (g, group) in enumerate(order_groups)
        partial = _iitm_repack_transition_matrix(𝐓s[g], group, nₘₐₓ)
        value_container .= value_container .+ partial
        for j in 1:pvars
            derivative_containers[j] .= derivative_containers[j] .+
                                        _iitm_repack_transition_matrix(∂𝐓s_by_var[j][g],
                group,
                nₘₐₓ)
        end
    end

    value = TransitionMatrix{CT, nₘₐₓ, typeof(value_container)}(value_container)
    jacobian = [TransitionMatrix{CT, nₘₐₓ, typeof(container)}(container)
                for container in derivative_containers]

    return LinearizationResult(value, jacobian, variables;
        metadata = (; backend = :iitm_nfold_analytic,
            variant = :nfold,
            λ = input.λ,
            nₘₐₓ,
            Nr,
            Nϑ = Nϑ_config,
            Nφ))
end

function supports_linearization(problem::LinearizationProblem,
        backend::IITMLinearization;
        output::Symbol = :transition_matrix,
        config = nothing)
    output == :transition_matrix ||
        return LinearizationSupport(false,
            "IITM analytical linearization only supports transition matrices")
    backend.variant in (:auto, :axisymmetric, :nfold, :arbitrary) ||
        return LinearizationSupport(false,
            "IITM analytical linearization variant must be :auto, :axisymmetric, :nfold, or :arbitrary")

    input = try
        _iitm_linearization_input(problem, config)
    catch err
        return LinearizationSupport(false, "failed to rebuild IITM input: $err")
    end
    isnothing(input) &&
        return LinearizationSupport(false,
            "IITM analytical linearization requires shape, λ, nₘₐₓ, Nr, and Nϑ")
    hasproperty(input.shape, :m) ||
        return LinearizationSupport(false,
            "IITM fixed-geometry analytical linearization requires shapes with an `m` material field")
    variant = _iitm_effective_variant(input, backend)
    if variant == :axisymmetric
        input.shape isa AbstractAxisymmetricShape ||
            return LinearizationSupport(false,
                "IITM axisymmetric analytical linearization requires an axisymmetric shape")
    elseif variant == :nfold
        input.shape isa AbstractNFoldShape ||
            return LinearizationSupport(false,
                "IITM nfold analytical linearization requires an n-fold shape")
        isnothing(input.Nφ) &&
            return LinearizationSupport(false,
                "IITM nfold analytical linearization requires Nφ")
    elseif variant == :arbitrary
        input.shape isa AbstractShape ||
            return LinearizationSupport(false,
                "IITM arbitrary analytical linearization requires an AbstractShape")
        isnothing(input.Nφ) &&
            return LinearizationSupport(false,
                "IITM arbitrary analytical linearization requires Nφ")
    else
        return LinearizationSupport(false,
            "IITM analytical linearization variant must be :auto, :axisymmetric, :nfold, or :arbitrary")
    end
    _linearization_variables_supported(variables(problem),
        _IITM_FIXED_GEOMETRY_LINEARIZATION_VARIABLES) ||
        return LinearizationSupport(false,
            "IITM fixed-geometry analytical linearization supports unique canonical variables drawn from $(_iitm_variable_list_message(_IITM_FIXED_GEOMETRY_LINEARIZATION_VARIABLES))")

    return LinearizationSupport(true, "")
end

function _checked_iitm_linearization_input(problem::LinearizationProblem,
        backend::IITMLinearization,
        config)
    support = supports_linearization(problem, backend; output = :transition_matrix,
        config)
    Bool(support) ||
        throw(UnsupportedLinearization(backend, :transition_matrix, support.reason))
    return _iitm_linearization_input(problem, config)
end

function linearize_transition_matrix(problem::LinearizationProblem,
        backend::IITMLinearization; config = nothing)
    input = _checked_iitm_linearization_input(problem, backend, config)
    variant = _iitm_effective_variant(input, backend)
    if variant == :axisymmetric
        return _iitm_axisymmetric_fixed_geometry_linearization(input, variables(problem))
    elseif variant == :nfold
        return _iitm_nfold_fixed_geometry_linearization(input, variables(problem))
    else
        return _iitm_arbitrary_fixed_geometry_linearization(input, variables(problem))
    end
end
