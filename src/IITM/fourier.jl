# Solve 𝐌 \ 𝐗 (= inv(𝐌) * 𝐗) through a factorization instead of forming the
# explicit inverse and an extra matrix product — cheaper and numerically
# cleaner on the IITM radial recursion's hot path. Arblib (Arb/Acb) matrices
# have no generic `\`/`lu`, so they fall back to their dedicated `inv`.
_iitm_ldiv(𝐌::AbstractMatrix{<:Union{Arb, Acb}}, 𝐗) = inv(𝐌) * 𝐗
_iitm_ldiv(𝐌, 𝐗) = 𝐌 \ 𝐗

# Capability predicate: true when the FFT path is available for this complex type.
# ComplexF64 → FFTW (hard dep); Acb → Arblib.dft! (hard dep). Other Complex{<:AbstractFloat}
# (Double64/BigFloat) become capable only when GenericFFT is loaded — see
# `ext/TransitionMatricesGenericFFTExt.jl`; without it they use the direct azimuthal sum.
_iitm_fft_capable(::Type{ComplexF64}) = true  # FFTW
_iitm_fft_capable(::Type{Acb}) = true         # Arblib.dft!
_iitm_fft_capable(::Type) = false             # direct fallback (incl. generic float w/o GenericFFT)

struct _AzimuthalFourierWorkspace{CT, P}
    contrast::Matrix{CT}
    contrast_inv::Matrix{CT}
    spectrum::Matrix{CT}
    spectrum_inv::Matrix{CT}
    plan::P
    coeff::OffsetArray{CT, 2, Matrix{CT}}
    coeff_inv::OffsetArray{CT, 2, Matrix{CT}}
end

# Acb workspace: holds Flint-native AcbVector/AcbMatrix scratch so that
# no Arb/Flint allocations happen inside the radial loop.
struct _AcbFourierWorkspace
    col_in::Arblib.AcbVector
    col_out::Arblib.AcbVector
    coeff_mat::AcbMatrix
    coeff_inv_mat::AcbMatrix
    coeff::OffsetArray{Acb, 2, AcbMatrix}
    coeff_inv::OffsetArray{Acb, 2, AcbMatrix}
    unit::Acb   # cached `1` at the working precision (reused each radial step)
end

# ---------- FFTW plan cache (ComplexF64 only) ----------
const _AZIMUTHAL_FFT_PLAN_CACHE = Dict{Tuple{Int, Int}, Any}()
const _AZIMUTHAL_FFT_PLAN_CACHE_LOCK = ReentrantLock()

function _azimuthal_fft_plan(::Type{ComplexF64}, Nφ, Nϑ)
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

# The GenericFFT plan for non-Float64 complex-float types lives in
# `ext/TransitionMatricesGenericFFTExt.jl` (GenericFFT is a weak dependency).

# ---------- Workspace constructors ----------

function _azimuthal_fourier_workspace(Nφ, Nϑ, nₘₐₓ; prec = 0)
    _azimuthal_fourier_workspace(ComplexF64, Nφ, Nϑ, nₘₐₓ)
end

function _azimuthal_fourier_workspace(::Type{CT}, Nφ, Nϑ, nₘₐₓ; prec = 0) where {CT <: Complex{<:AbstractFloat}}
    contrast = zeros(CT, Nφ, Nϑ)
    contrast_inv = zeros(CT, Nφ, Nϑ)
    spectrum = zeros(CT, Nφ, Nϑ)
    spectrum_inv = zeros(CT, Nφ, Nϑ)
    plan = _azimuthal_fft_plan(CT, Nφ, Nϑ)
    coeff = OffsetArray(zeros(CT, 4nₘₐₓ + 1, Nϑ), (-2nₘₐₓ):(2nₘₐₓ), 1:Nϑ)
    coeff_inv = OffsetArray(zeros(CT, 4nₘₐₓ + 1, Nϑ), (-2nₘₐₓ):(2nₘₐₓ), 1:Nϑ)
    return _AzimuthalFourierWorkspace{CT, typeof(plan)}(
        contrast, contrast_inv, spectrum, spectrum_inv, plan, coeff, coeff_inv)
end

# Acb workspace: preallocated AcbVector/AcbMatrix scratch; no Flint allocations
# inside the radial loop.
function _azimuthal_fourier_workspace(::Type{Acb}, Nφ, Nϑ, nₘₐₓ; prec)
    col_in = Arblib.AcbVector(Nφ; prec = prec)
    col_out = Arblib.AcbVector(Nφ; prec = prec)
    coeff_mat = AcbMatrix(4nₘₐₓ + 1, Nϑ; prec = prec)
    coeff_inv_mat = AcbMatrix(4nₘₐₓ + 1, Nϑ; prec = prec)
    coeff = OffsetArray(coeff_mat, (-2nₘₐₓ):(2nₘₐₓ), 1:Nϑ)
    coeff_inv = OffsetArray(coeff_inv_mat, (-2nₘₐₓ):(2nₘₐₓ), 1:Nϑ)
    unit = Acb(1; prec = prec)
    return _AcbFourierWorkspace(col_in, col_out, coeff_mat, coeff_inv_mat, coeff, coeff_inv, unit)
end

# ---------- Mode-bin helpers ----------

function _azimuthal_fourier_mode_bins(nₘₐₓ, period = 1)
    qs = [q for q in (-2nₘₐₓ):(2nₘₐₓ) if q % period == 0]
    bins = [q ÷ period for q in qs]
    return qs, bins
end

# ---------- Forward column DFT: ComplexF64 (FFTW) ----------
# FFTW plans support mul!(out, plan, in) in-place via AbstractFFTs.
function _apply_forward_dft!(spectrum::Matrix{ComplexF64},
        contrast::Matrix{ComplexF64},
        plan)
    mul!(spectrum, plan, contrast)
end

# The generic (non-ComplexF64) forward column DFT — used by the GenericFFT plan — lives in
# `ext/TransitionMatricesGenericFFTExt.jl` (GenericFFT is a weak dependency).

# ---------- Fourier coefficient extraction: generic Complex{<:AbstractFloat} ----------

function _azimuthal_fourier_coefficients(ε::AbstractMatrix{CT}, nₘₐₓ, wφ,
        workspace::_AzimuthalFourierWorkspace{CT},
        mode_bins) where {CT <: Complex{<:AbstractFloat}}
    Nφ, Nϑ = size(ε)
    qs, bins = mode_bins

    @. workspace.contrast = ε - 1
    @. workspace.contrast_inv = workspace.contrast / ε
    _apply_forward_dft!(workspace.spectrum, workspace.contrast, workspace.plan)
    _apply_forward_dft!(workspace.spectrum_inv, workspace.contrast_inv, workspace.plan)

    # Write into preallocated workspace arrays (no allocation).
    # Alias safety: callers consume coeff/coeff_inv before the next radial step.
    coeff = workspace.coeff
    coeff_inv = workspace.coeff_inv

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

# ---------- Fourier coefficient extraction: Acb (Arblib.dft! backend) ----------
# Uses AcbMatrix (Arblib native) for the coeff arrays — Julia's Matrix{Acb} is
# unsafe due to GC interaction with the Flint memory model.
# Precision is inferred from the ε matrix elements at runtime.

function _azimuthal_fourier_coefficients(ε::AbstractMatrix{Acb}, nₘₐₓ, wφ,
        workspace::_AcbFourierWorkspace,
        mode_bins)
    Nφ, Nϑ = size(ε)
    qs, bins = mode_bins
    prec = Arblib.precision(ε[1, 1])

    # Reuse preallocated AcbVector/AcbMatrix scratch + cached unit (no Flint allocations).
    # Alias safety: callers consume coeff/coeff_inv before the next radial step.
    col_in = workspace.col_in
    col_out = workspace.col_out
    coeff = workspace.coeff
    coeff_inv = workspace.coeff_inv
    one_prec = workspace.unit

    for i in 1:Nϑ
        # contrast column: ε[j,i] - 1
        for j in 1:Nφ
            col_in[j] = ε[j, i] - one_prec
        end
        Arblib.dft!(col_out, col_in; len = Nφ, prec = prec)
        for iq in eachindex(qs)
            q = qs[iq]
            idx = mod(-bins[iq], Nφ) + 1
            coeff[q, i] = wφ * col_out[idx]
        end

        # contrast_inv column: (ε[j,i] - 1) / ε[j,i]
        for j in 1:Nφ
            col_in[j] = (ε[j, i] - one_prec) / ε[j, i]
        end
        Arblib.dft!(col_out, col_in; len = Nφ, prec = prec)
        for iq in eachindex(qs)
            q = qs[iq]
            idx = mod(-bins[iq], Nφ) + 1
            coeff_inv[q, i] = wφ * col_out[idx]
        end
    end

    return coeff, coeff_inv
end
