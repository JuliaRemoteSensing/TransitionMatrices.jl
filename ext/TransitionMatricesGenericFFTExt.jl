module TransitionMatricesGenericFFTExt

# Generic-type azimuthal FFT for the n-fold IITM, loaded on demand. GenericFFT is a *weak*
# dependency: ComplexF64 uses FFTW and Acb uses Arblib.dft! (both hard deps), so only the
# Complex{Double64} / Complex{BigFloat} FFT path needs GenericFFT. Without it, those types
# fall back to the direct azimuthal sum (slower but correct) — see `_iitm_fft_capable`.

using TransitionMatrices
using GenericFFT: plan_fft

# Make non-Float64 complex-float types FFT-capable. ComplexF64 keeps its more-specific core
# method (FFTW); Acb keeps Arblib.dft!.
TransitionMatrices._iitm_fft_capable(::Type{<:Complex{<:AbstractFloat}}) = true

# Fresh GenericFFT plan each call: a BigFloat plan bakes in precision-specific twiddles and
# must not be cached/reused across precision changes (no FFTW flags here).
function TransitionMatrices._azimuthal_fft_plan(::Type{CT}, Nφ, Nϑ) where {CT <: Complex{<:AbstractFloat}}
    return plan_fft(zeros(CT, Nφ, Nϑ), 1)
end

# GenericFFT's plan supports `plan * A` (allocating), not in-place mul!.
function TransitionMatrices._apply_forward_dft!(spectrum::Matrix{CT}, contrast::Matrix{CT},
        plan) where {CT <: Complex{<:AbstractFloat}}
    spectrum .= plan * contrast
end

end
