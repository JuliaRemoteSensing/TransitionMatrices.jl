struct _AzimuthalFourierWorkspace{P}
    contrast::Matrix{ComplexF64}
    contrast_inv::Matrix{ComplexF64}
    spectrum::Matrix{ComplexF64}
    spectrum_inv::Matrix{ComplexF64}
    plan::P
end

const _AZIMUTHAL_FFT_PLAN_CACHE = Dict{Tuple{Int, Int}, Any}()
const _AZIMUTHAL_FFT_PLAN_CACHE_LOCK = ReentrantLock()

function _azimuthal_fft_plan(Nφ, Nϑ)
    key = (Nφ, Nϑ)
    lock(_AZIMUTHAL_FFT_PLAN_CACHE_LOCK)
    try
        return get!(_AZIMUTHAL_FFT_PLAN_CACHE, key) do
            FFTW.plan_fft(zeros(ComplexF64, Nφ, Nϑ), 1;
                          flags = FFTW.ESTIMATE | FFTW.UNALIGNED)
        end
    finally
        unlock(_AZIMUTHAL_FFT_PLAN_CACHE_LOCK)
    end
end

function _azimuthal_fourier_workspace(Nφ, Nϑ)
    contrast = Matrix{ComplexF64}(undef, Nφ, Nϑ)
    contrast_inv = similar(contrast)
    spectrum = similar(contrast)
    spectrum_inv = similar(contrast)
    plan = _azimuthal_fft_plan(Nφ, Nϑ)
    return _AzimuthalFourierWorkspace(contrast, contrast_inv, spectrum, spectrum_inv,
                                      plan)
end

function _azimuthal_fourier_mode_bins(nₘₐₓ, period = 1)
    qs = [q for q in (-2nₘₐₓ):(2nₘₐₓ) if q % period == 0]
    bins = [q ÷ period for q in qs]
    return qs, bins
end

function _azimuthal_fourier_coefficients(ε::AbstractMatrix{ComplexF64}, nₘₐₓ, wφ,
                                         workspace::_AzimuthalFourierWorkspace,
                                         mode_bins)
    Nφ, Nϑ = size(ε)
    qs, bins = mode_bins

    @. workspace.contrast = ε - 1
    @. workspace.contrast_inv = workspace.contrast / ε
    mul!(workspace.spectrum, workspace.plan, workspace.contrast)
    mul!(workspace.spectrum_inv, workspace.plan, workspace.contrast_inv)

    coeff = OffsetArray(zeros(ComplexF64, 4nₘₐₓ + 1, Nϑ), (-2nₘₐₓ):(2nₘₐₓ),
                        1:Nϑ)
    coeff_inv = similar(coeff)

    for i in 1:Nϑ
        for iq in eachindex(qs)
            q = qs[iq]
            idx = mod(-bins[iq], Nφ) + 1
            coeff[q, i] = wφ * workspace.spectrum[idx, i]
            coeff_inv[q, i] = wφ * workspace.spectrum_inv[idx, i]
        end
    end

    return coeff, coeff_inv
end
