function routine_mishchenko(threshold, ndgs, nₘₐₓ_only)
    Qsca0 = 0.0
    Qext0 = 0.0
    nₘₐₓ_converged = false

    function routine(nₘₐₓ, Ng, Qsca, Qext)
        Qsca = Float64(Qsca)
        Qext = Float64(Qext)
        ΔQsca = abs((Qsca0 - Qsca) / Qsca)
        ΔQext = abs((Qext0 - Qext) / Qext)
        Qsca0 = Qsca
        Qext0 = Qext
        Δ = max(ΔQsca, ΔQext)

        if !nₘₐₓ_converged
            @debug "nₘₐₓ iteration: nₘₐₓ=$nₘₐₓ Ng=$Ng ΔQsca=$ΔQsca ΔQext=$ΔQext"
            if Δ < threshold
                @debug "nₘₐₓ converged"
                nₘₐₓ_converged = true

                if nₘₐₓ_only
                    return -1, -1
                end
                return nₘₐₓ, Ng + 4
            else
                if nₘₐₓ_only
                    return nₘₐₓ + 1, Ng
                end
                return nₘₐₓ + 1, Ng + ndgs
            end
        else
            @debug "Ng iteration: nₘₐₓ=$nₘₐₓ Ng=$Ng ΔQsca=$ΔQsca ΔQext=$ΔQext"
            if Δ < threshold
                @debug "Ng converged"
                return -1, -1
            else
                return nₘₐₓ, Ng + 4
            end
        end
    end

    return routine
end

@doc raw"""
```
transition_matrix(s::AbstractAxisymmetricShape{T, CT}, λ; ...) where {T, CT}
```

Main function for calculating the transition matrix of an axisymmetric shape.

Parameters:

- `s`: Axisymmetric shape
- `λ`: Wavelength

Keyword arguments:

- `threshold`: Convergence threshold, defaults to `0.0001`
- `ndgs`: Number of discrete Gauss quadrature points per degree, defaults to `4`
- `routine_generator`: Function that generates the convergence routine, defaults to `routine_mishchenko`
- `nₛₜₐᵣₜ`: The initial number of the truncation order, defaults to `0`, and the actual value will be calculated automatically using the following formula:

```math
n_{\text{start}} = \max(4, \lceil kr_{\max} + 4.05 \sqrt[3]{kr_{\max}} \rceil)
```

- `Ngₛₜₐᵣₜ`: The initial number of discrete Gauss quadrature points, defaults to `nₛₜₐᵣₜ * ndgs`
- `nₘₐₓ_only`: Only check convergence of `nₘₐₓ`, defaults to `false`, meaning that convergence of both `nₘₐₓ` and `Ng` will be checked
- `full`: Whether to calculate the whole transition matrix in each iteration, defaults to `false`, and only the `m=m′=0` block will be calculated
- `reuse`: Whether to reuse the cache of the previous iteration, defaults to `true`
- `maxiter`: Maximum number of iterations, defaults to `20`
- `zerofn`: Function that defines the type used for numerical integration, defaults to `() -> zero(CT)` where `CT` is defined by the shape
"""
function transition_matrix(s::AbstractAxisymmetricShape{T, CT}, λ; threshold = 0.0001,
        ndgs = 4, routine_generator = routine_mishchenko,
        nₛₜₐᵣₜ = 0, Ngₛₜₐᵣₜ = nₛₜₐᵣₜ * ndgs, nₘₐₓ_only = false,
        full = false, reuse = true, maxiter = 20,
        zerofn = () -> zero(CT)) where {T, CT}
    if nₛₜₐᵣₜ == 0
        kr = 2π * rmax(s) / λ
        if nₛₜₐᵣₜ == 0
            nₛₜₐᵣₜ = max(4, ceil(Int, kr + 4.05 * ∛kr))
            Ngₛₜₐᵣₜ = nₛₜₐᵣₜ * ndgs
        end
    end

    routine = routine_generator(threshold, ndgs, nₘₐₓ_only)
    nₘₐₓ, Ng = nₛₜₐᵣₜ, Ngₛₜₐᵣₜ
    counter = 0
    while true
        cache = nothing
        counter += 1

        if !full
            T₀, cache = transition_matrix_m₀(s, λ, nₘₐₓ, Ng; reuse = true)
            Qext, Qsca = extinction_efficiency_m₀(T₀), scattering_efficiency_m₀(T₀)
        else
            𝐓 = transition_matrix(s, λ, nₘₐₓ, Ng)
            Qext, Qsca = extinction_cross_section(𝐓, λ), scattering_cross_section(𝐓, λ)
        end

        @debug "Qsca = $Qsca, Qext = $Qext"
        nₘₐₓ′, Ng′ = routine(nₘₐₓ, Ng, Qsca, Qext)
        if nₘₐₓ′ == -1
            if full
                return 𝐓
            else
                Ts = [T₀]
                for m in 1:nₘₐₓ
                    @debug "Calculating for m = $m"
                    Tₘ = transition_matrix_m(m, s, λ, nₘₐₓ, Ng;
                        cache = reuse ? cache : nothing)
                    push!(Ts, Tₘ)
                end
                return AxisymmetricTransitionMatrix{CT, nₘₐₓ, typeof(Ts), T}(Ts)
            end
        else
            nₘₐₓ, Ng = nₘₐₓ′, Ng′
        end

        if counter > maxiter
            error("Maximum number of iterations reached, failed to converge.")
        end
    end
end
