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
