@doc raw"""
A general T-Matrix ``T_{m n m^{\prime} n^{\prime}}^{k l}`` stored in a 6-dimensional array, in the order ``(m, n, m^{\prime}, n^{\prime}, k, l)``.
"""
abstract type AbstractTransitionMatrix{CT, N} <: AbstractArray{CT, 6} end

"""
Get the maximum order of a T-Matrix.
"""
order(::AbstractTransitionMatrix{CT, N}) where {CT, N} = N

Base.size(::AbstractTransitionMatrix{CT, N}) where {CT, N} = (2N + 1, N, 2N + 1, N, 2, 2)
function Base.axes(::AbstractTransitionMatrix{CT, N}) where {CT, N}
    ((-N):N, 1:N, (-N):N, 1:N, 1:2, 1:2)
end

"""
Concrete type for a general T-Matrix.
"""
struct TransitionMatrix{CT, N, V <: AbstractArray{CT}} <: AbstractTransitionMatrix{CT, N}
    container::V
end

Base.getindex(tm::TransitionMatrix{CT, N}, idx) where {CT, N} = getindex(tm.container, idx)
function Base.getindex(tm::TransitionMatrix{CT, N}, idxs...) where {CT, N}
    getindex(tm.container, idxs...)
end

@doc raw"""
```
rotate(𝐓::AbstractTransitionMatrix{CT, N}, rot::Rotation{3})
```

Rotate the given T-Matrix `𝐓` by the Euler angle `rot` and generate a new T-Matrix.

For a general T-Matrix, Eq. (5.29) in Mishchenko et al. (2002) is used as a fallback. A `TransitionMatrix` will be returned, which is the most general yet concrete type.

```math
T_{m n m^{\prime} n^{\prime}}^{p p′}(L ; \alpha, \beta, \gamma)=\sum_{m_1=-n}^n \sum_{m_2=-n^{\prime}}^{n^{\prime}} D_{m m_1}^n(\alpha, \beta, \gamma) T_{m_1 n m_2 n^{\prime}}^{p p′}(P) D_{m_2 m^{\prime}}^{n^{\prime}}(-\gamma,-\beta,-\alpha)\quad p,p′=1,2
```
"""
function rotate(𝐓::AbstractTransitionMatrix{CT, N}, rot::Rotation{3}) where {CT, N}
    # Get the Euler angle in Z-Y-Z order.
    zyz = RotZYZ(rot)
    α, β, γ = zyz.theta1, zyz.theta2, zyz.theta3

    # Calculate the wigner-d functions that will be used.
    d = OffsetArray(zeros(CT, 2N + 1, 2N + 1, N + 1), (-N):N, (-N):N, 0:N)
    for m in (-N):N
        for m′ in (-N):N
            sₘᵢₙ = max(abs(m), abs(m′))
            wigner_d_recursion!(view(d, m, m′, sₘᵢₙ:N), m, m′, N, β)
        end
    end

    # Calculate the coefficients used for wigner-D functions
    coeff = OffsetArray([cis(-(m * α + m′ * γ)) for m in (-N):N, m′ in (-N):N], (-N):N,
                        (-N):N)

    # Calculate the rotated T-Matrix
    𝐓′ = similar(𝐓)
    fill!(𝐓′, 0)

    # Enable multi-threading
    Threads.@threads for (n′, m′) in OrderDegreeIterator(N)
        for p in 1:2, p′ in 1:2
            for (n, m) in OrderDegreeIterator(N)
                for m₂ in (-n′):n′, m₁ in (-n):n
                    sig = iseven(m′ + m₂) ? 1 : -1
                    𝐓′[m, n, m′, n′, p, p′] += coeff[m, m₁] * d[m, m₁, n] *
                                               conj(coeff[m′, m₂]) * d[m₂, m′, n′] * sig *
                                               𝐓[m₁, n, m₂, n′, p, p′]
                end
            end
        end
    end

    TransitionMatrix{CT, N, typeof(𝐓′)}(𝐓′)
end

@doc raw"""
```
amplitude_matrix(𝐓::AbstractTransitionMatrix{CT, N}, ϑᵢ, φᵢ, ϑₛ, φₛ; λ=2π)
```

Calculate the amplitude matrix of the given T-Matrix `𝐓` at the given incidence and scattering angles. 

Parameters:

- `𝐓`: the T-Matrix of the scatterer.
- `ϑᵢ`: the incidence zenith angle in radians.
- `φᵢ`: the incidence azimuthal angle in radians.
- `ϑₛ`: the scattering zenith angle in radians.
- `φₛ`: the scattering azimuthal angle in radians.
- `λ`: the wavelength of the incident wave in the host medium. Default to 2π.

For a general T-Matrix, Eq. (5.11) -- Eq. (5.17) in Mishchenko et al. (2002) is used as a fallback.

```math
\begin{array}{l}
S_{11}\left(\hat{\mathbf{n}}^{\text {sca }}, \hat{\mathbf{n}}^{\text {inc }}\right)=\frac{1}{k_1} \sum_{n=1}^{\infty} \sum_{n^{\prime}=1}^{\infty} \sum_{m=-n}^n \sum_{m^{\prime}=-n^{\prime}}^{n^{\prime}} \alpha_{m n m^{\prime} n^{\prime}}\left[T_{m n m^{\prime} n^{\prime}}^{11} \pi_{m n}\left(\vartheta^{\text {sca }}\right) \pi_{m^{\prime} n^{\prime}}\left(\vartheta^{\text {inc }}\right)\right. \\
+T_{m n m^{\prime} n^{\prime}}^{21} \tau_{m n}\left(\vartheta^{\text {sca }}\right) \pi_{m^{\prime} n^{\prime}}\left(\vartheta^{\mathrm{inc}}\right)+T_{m n m^{\prime} n^{\prime}}^{12} \pi_{m n}\left(\vartheta^{\text {sca }}\right) \tau_{m^{\prime} n^{\prime}}\left(\vartheta^{\mathrm{inc}}\right) \\
\left.+T_{m n m^{\prime} n^2}^{22} \tau_{m n}\left(\vartheta^{\text {sca }}\right) \tau_{m^{\prime} n^{\prime}}\left(\vartheta^{\text {inc }}\right)\right] \exp \left[\mathrm{i}\left(m \varphi^{\text {sca }}-m^{\prime} \varphi^{\text {inc }}\right)\right] \text {, } \\
S_{12}\left(\hat{\mathbf{n}}^{\mathrm{sca}}, \hat{\mathbf{n}}^{\mathrm{inc}}\right)=\frac{1}{\mathrm{i} k_1} \sum_{n=1}^{\infty} \sum_{n^{\prime}=1}^{\infty} \sum_{m=-n}^n \sum_{m^{\prime}=-n^{\prime}}^{n^{\prime}} \alpha_{m n m^{\prime} n^{\prime}}\left[T_{m n m^{\prime} n^{\prime}}^{11} \pi_{m n}\left(\vartheta^{\mathrm{sca}}\right) \tau_{m^{\prime} n^{\prime}}\left(\vartheta^{\mathrm{inc}}\right)\right. \\
+T_{m n m^{\prime} n^{\prime}}^{21} \tau_{m n}\left(\vartheta^{\text {sca }}\right) \tau_{m^{\prime} n^{\prime}}\left(\vartheta^{\mathrm{inc}}\right)+T_{m n m^{\prime} n^{\prime}}^{12} \pi_{m n}\left(\vartheta^{\text {sca }}\right) \pi_{m^{\prime} n^{\prime}}\left(\vartheta^{\text {inc }}\right) \\
\left.+T_{m n m^{\prime} n^{\prime}}^{22} \tau_{m n}\left(\vartheta^{\text {sca }}\right) \pi_{m^{\prime} n^{\prime}}\left(\vartheta^{\text {inc }}\right)\right] \exp \left[\mathrm{i}\left(m \varphi^{\text {sca }}-m^{\prime} \varphi^{\text {inc }}\right)\right] \text {, } \\
S_{21}\left(\hat{\mathbf{n}}^{\mathrm{sca}}, \hat{\mathbf{n}}^{\mathrm{inc}}\right)=\frac{\mathrm{i}}{k_1} \sum_{n=1}^{\infty} \sum_{n^{\prime}=1}^{\infty} \sum_{m=-n}^n \sum_{m^{\prime}=-n^{\prime}}^{n^{\prime}} \alpha_{m n m^{\prime} n^{\prime}}\left[T_{m n m^{\prime} n^{\prime}}^{11} \tau_{m n}\left(\vartheta^{\mathrm{sca}}\right) \pi_{m^{\prime} n^{\prime}}\left(\vartheta^{\mathrm{inc}}\right)\right. \\
+T_{m n m^{\prime} n^{\prime}}^{21} \pi_{m n}\left(\vartheta^{\text {sca }}\right) \pi_{m^{\prime} n^{\prime}}\left(\vartheta^{\mathrm{inc}}\right)+T_{m n m^{\prime} n^{\prime}}^{12} \tau_{m n}\left(\vartheta^{\text {sca }}\right) \tau_{m^{\prime} n^{\prime}}\left(\vartheta^{\text {inc }}\right) \\
\left.+T_{m n m^{\prime} n^{\prime}}^{22} \pi_{m n}\left(\vartheta^{\text {sca }}\right) \tau_{m^{\prime} n^{\prime}}\left(\vartheta^{\text {inc }}\right)\right] \exp \left[\mathrm{i}\left(m \varphi^{\text {sca }}-m^{\prime} \varphi^{\text {inc }}\right)\right] \text {, } \\
S_{22}\left(\hat{\mathbf{n}}^{\text {sca }}, \hat{\mathbf{n}}^{\mathrm{inc}}\right)=\frac{1}{k_1} \sum_{n=1}^{\infty} \sum_{n^{\prime}=1}^{\infty} \sum_{m=-n}^n \sum_{m^{\prime}=-n^{\prime}}^{n^{\prime}} \alpha_{m n m^{\prime} n^{\prime}}\left[T_{m n n^{\prime} n^{\prime}}^{11} \tau_{m n}\left(\vartheta^{\text {sca }}\right) \tau_{m^{\prime} n^{\prime}}\left(\vartheta^{\text {inc }}\right)\right. \\
+T_{m n m^{\prime} n^{\prime}}^{21} \pi_{m n}\left(\vartheta^{\text {sca }}\right) \tau_{m^{\prime} n^{\prime}}\left(\vartheta^{\text {inc }}\right)+T_{m n m^{\prime} n^{12}}^{12} \tau_{m n}\left(\vartheta^{\text {sca }}\right) \pi_{m^{\prime} n^{\prime}}\left(\vartheta^{\text {inc }}\right) \\
\left.+T_{m n m^{\prime} n^{\prime}}^{22} \pi_{m n}\left(\vartheta^{\text {sca }}\right) \pi_{m^{\prime} n^{\prime}}\left(\vartheta^{\text {inc }}\right)\right] \exp \left[\mathrm{i}\left(m \varphi^{\text {sca }}-m^{\prime} \varphi^{\text {inc }}\right)\right], \\
\end{array}
```

Where

```math
\begin{array}{l}
\alpha_{m n m^{\prime} n^{\prime}}=\mathrm{i}^{n^{\prime}-n-1}(-1)^{m+m^{\prime}}\left[\frac{(2 n+1)\left(2 n^{\prime}+1\right)}{n(n+1) n^{\prime}\left(n^{\prime}+1\right)}\right]^{1 / 2}, \\
\pi_{m n}(\vartheta)=\frac{m d_{0 m}^n(\vartheta)}{\sin \vartheta}, \quad \pi_{-m n}(\vartheta)=(-1)^{m+1} \pi_{m n}(\vartheta), \\
\tau_{m n}(\vartheta)=\frac{\mathrm{d} d_{0 m}^n(\vartheta)}{\mathrm{d} \vartheta}, \quad \tau_{-m n}(\vartheta)=(-1)^m \tau_{m n}(\vartheta)
\end{array}
```
"""
function amplitude_matrix(𝐓::AbstractTransitionMatrix{CT, N}, ϑᵢ, φᵢ, ϑₛ, φₛ;
                          λ = 2π) where {CT, N}
    T = real(CT)
    k₁ = 2π / λ
    𝐒₁₁, 𝐒₁₂, 𝐒₂₁, 𝐒₂₂ = zero(CT), zero(CT), zero(CT), zero(CT)

    πᵢ = OffsetArray(zeros(T, 2N + 1, N + 1), (-N):N, 0:N)
    τᵢ = OffsetArray(zeros(T, 2N + 1, N + 1), (-N):N, 0:N)
    πₛ = OffsetArray(zeros(T, 2N + 1, N + 1), (-N):N, 0:N)
    τₛ = OffsetArray(zeros(T, 2N + 1, N + 1), (-N):N, 0:N)
    for m in 0:N
        wigner_d_recursion!(view(πᵢ, m, m:N),
                            0, m, N, ϑᵢ;
                            deriv = view(τᵢ, m, m:N))

        wigner_d_recursion!(view(πₛ, m, m:N),
                            0, m, N, ϑₛ;
                            deriv = view(τₛ, m, m:N))
    end

    for n in 1:N
        for m in 0:n
            πᵢ[m, n] = pi_func(T, m, n, ϑᵢ; d = πᵢ[m, n])
            πₛ[m, n] = pi_func(T, m, n, ϑₛ; d = πₛ[m, n])
            if m > 0
                πᵢ[-m, n] = (-1)^((m + 1) & 1) * πᵢ[m, n]
                πₛ[-m, n] = (-1)^((m + 1) & 1) * πₛ[m, n]
                τᵢ[-m, n] = (-1)^(m & 1) * τᵢ[m, n]
                τₛ[-m, n] = (-1)^(m & 1) * τₛ[m, n]
            end
        end
    end

    for n′ in 1:N, n in 1:N
        αₙ = 1.0im^((n′ - n - 1) & 3) *
             √(T(2n + 1) * (2n′ + 1) / (n * (n + 1) * n′ * (n′ + 1)))
        for m′ in (-n′):n′
            for m in (-n):n
                α = (-1.0)^((m + m′) & 1) * αₙ
                expiφ = cis(m * φₛ - m′ * φᵢ)
                𝐒₁₁ += α *
                       (𝐓[m, n, m′, n′, 1, 1] * πₛ[m, n] * πᵢ[m′, n′] +
                        𝐓[m, n, m′, n′, 1, 2] * πₛ[m, n] * τᵢ[m′, n′] +
                        𝐓[m, n, m′, n′, 2, 1] * τₛ[m, n] * πᵢ[m′, n′] +
                        𝐓[m, n, m′, n′, 2, 2] * τₛ[m, n] * τᵢ[m′, n′]) * expiφ

                𝐒₁₂ += α *
                       (𝐓[m, n, m′, n′, 1, 1] * πₛ[m, n] * τᵢ[m′, n′] +
                        𝐓[m, n, m′, n′, 1, 2] * πₛ[m, n] * πᵢ[m′, n′] +
                        𝐓[m, n, m′, n′, 2, 1] * τₛ[m, n] * τᵢ[m′, n′] +
                        𝐓[m, n, m′, n′, 2, 2] * τₛ[m, n] * πᵢ[m′, n′]) * expiφ

                𝐒₂₁ += α *
                       (𝐓[m, n, m′, n′, 1, 1] * τₛ[m, n] * πᵢ[m′, n′] +
                        𝐓[m, n, m′, n′, 1, 2] * τₛ[m, n] * τᵢ[m′, n′] +
                        𝐓[m, n, m′, n′, 2, 1] * πₛ[m, n] * πᵢ[m′, n′] +
                        𝐓[m, n, m′, n′, 2, 2] * πₛ[m, n] * τᵢ[m′, n′]) * expiφ

                𝐒₂₂ += α *
                       (𝐓[m, n, m′, n′, 1, 1] * τₛ[m, n] * τᵢ[m′, n′] +
                        𝐓[m, n, m′, n′, 1, 2] * τₛ[m, n] * πᵢ[m′, n′] +
                        𝐓[m, n, m′, n′, 2, 1] * πₛ[m, n] * τᵢ[m′, n′] +
                        𝐓[m, n, m′, n′, 2, 2] * πₛ[m, n] * πᵢ[m′, n′]) * expiφ
            end
        end
    end

    return (@SMatrix [𝐒₁₁ 𝐒₁₂/1im; 𝐒₂₁*1im 𝐒₂₂]) ./ k₁
end

@doc raw"""
```
orientation_average(𝐓::AbstractTransitionMatrix{CT, N}, pₒ; Nα = 10, Nβ = 10, Nγ = 10) where {CT, N}
```

Calculate the orientation average of a transition matrix using numerical integration, given the orientation distribution function ``p_o(\alpha,\beta,\gamma)``. 

```math
\langle T_{m n m^{\prime} n^{\prime}}^{p p^{\prime}}(L)\rangle = \int_0^{2\pi}\mathrm{d}\alpha\int_0^{\pi}\mathrm{d}\beta\sin\beta\int_0^{2\pi}\mathrm{d}\gamma p_o(\alpha,\beta,\gamma) T_{m n m^{\prime} n^{\prime}}^{p p^{\prime}}(L; \alpha,\beta,\gamma)
```

Parameters:

- `𝐓`: the T-Matrix to be orientation averaged.
- `pₒ`: the orientation distribution function. Note that the ``\sin\beta`` part is already included.
- `Nα`: the number of points used in the numerical integration of ``\alpha``. Default to 10.
- `Nβ`: the number of points used in the numerical integration of ``\beta``. Default to 10.
- `Nγ`: the number of points used in the numerical integration of ``\gamma``. Default to 10.

!!! note

    This is the fallback method and does not utilize any symmetry, so it is expected to be slow. You should use specified versions of this function, or implement your own if there is no suited version for your combination of T-Matrix and orientation distribution function.

    You may also need to test the convergence of `Nα`, `Nβ` and `Nγ` manually. If any one is too small, there will be large errors in the results.
"""
function orientation_average(𝐓::AbstractTransitionMatrix{CT, N}, pₒ; Nα = 10, Nβ = 10,
                             Nγ = 10) where {CT, N}
    T̄ = similar(𝐓)
    fill!(T̄, zero(CT))

    # Transform to [0, 2pi]
    xa, wa = gausslegendre(Nα)
    @. xa = (xa + 1) * π
    @. wa *= π

    # Integrate cos(beta) from -1 to 1
    xb, wb = gausslegendre(Nβ)
    # Get beta from cos(beta)
    xb = acos.(xb)

    # Transform to [0, 2pi]
    xc, wc = gausslegendre(Nγ)
    @. xc = (xc + 1) * π
    @. wc *= π

    for (α, wi) in zip(xa, wa)
        for (β, wj) in zip(xb, wb)
            for (γ, wk) in zip(xc, wc)
                𝐓αβγ = rotate(𝐓, RotZYZ(α, β, γ))
                @. T̄ += pₒ(α, β, γ) * 𝐓αβγ * wi * wj * wk
            end
        end
    end

    TransitionMatrix{CT, N, typeof(T̄)}(T̄)
end

@doc raw"""
```
scattering_cross_section(𝐓::AbstractTransitionMatrix{CT, N}, λ=2π) where {CT, N}
```

Calculate the scattering cross section per particle averaged over the uniform orientation distribution, according to Eq. (5.140) in Mishchenko et al. (2002).

```math
\left\langle C_{\mathrm{sca}}\right\rangle=\frac{2 \pi}{k_1^2} \sum_{n=1}^{\infty} \sum_{m=-n}^n \sum_{n^{\prime}=1}^{\infty} \sum_{m^{\prime}=-n^{\prime}}^{n^{\prime}} \sum_{k=1}^2 \sum_{l=1}^2\left|T_{m n m^{\prime} n^{\prime}}^{k l}(P)\right|^2
```

Parameters:

- `𝐓`: the T-Matrix of the scatterer.
- `λ`: the wavelength of the incident wave in the host medium. Default to 2π.
"""
function scattering_cross_section(𝐓::AbstractTransitionMatrix{CT, N},
                                  λ = 2π) where {CT, N}
    return sum(abs2, 𝐓) * λ^2 / 2π
end

@doc raw"""
```
extinction_cross_section(𝐓::AbstractTransitionMatrix{CT, N}, λ=2π) where {CT, N}
```

Calculate the extinction cross section per particle averaged over the uniform orientation distribution, according to Eq. (5.102) in Mishchenko et al. (2002).

```math
\left\langle C_{\mathrm{ext}}\right\rangle=-\frac{2 \pi}{k_1^2} \operatorname{Re} \sum_{n=1}^{\infty} \sum_{m=-n}^n\left[T_{m n n n}^{11}(P)+T_{m n m n}^{22}(P)\right]
```

Parameters:

- `𝐓`: the T-Matrix of the scatterer.
- `λ`: the wavelength of the incident wave in the host medium. Default to 2π.
"""
function extinction_cross_section(𝐓::AbstractTransitionMatrix{CT, N},
                                  λ = 2π) where {CT, N}
    Cᵉˣᵗ = zero(CT)

    for n in 1:N
        for m in (-n):n
            Cᵉˣᵗ += 𝐓[m, n, m, n, 1, 1] + 𝐓[m, n, m, n, 2, 2]
        end
    end

    -real(Cᵉˣᵗ) * λ^2 / 2π
end

"""
```
phase_matrix(𝐒::AbstractMatrix)
```

Calculate the phase matrix `𝐙` from the amplitude matrix `𝐒`, according to Eq. (2.106) -- Eq. (2.121) in Mishchenko et al. (2002).
"""
function phase_matrix(𝐒::AbstractMatrix)
    𝐙₁₁ = 0.5 * (𝐒[1, 1] * 𝐒[1, 1]' + 𝐒[1, 2] * 𝐒[1, 2]' + 𝐒[2, 1] * 𝐒[2, 1]' +
           𝐒[2, 2] * 𝐒[2, 2]')
    𝐙₁₂ = 0.5 * (𝐒[1, 1] * 𝐒[1, 1]' - 𝐒[1, 2] * 𝐒[1, 2]' + 𝐒[2, 1] * 𝐒[2, 1]' -
           𝐒[2, 2] * 𝐒[2, 2]')
    𝐙₁₃ = -𝐒[1, 1] * 𝐒[1, 2]' - 𝐒[2, 2] * 𝐒[2, 1]'
    𝐙₁₄ = 1.0im * (𝐒[1, 1] * 𝐒[1, 2]' - 𝐒[2, 2] * 𝐒[2, 1]')
    𝐙₂₁ = 0.5 * (𝐒[1, 1] * 𝐒[1, 1]' + 𝐒[1, 2] * 𝐒[1, 2]' - 𝐒[2, 1] * 𝐒[2, 1]' -
           𝐒[2, 2] * 𝐒[2, 2]')
    𝐙₂₂ = 0.5 * (𝐒[1, 1] * 𝐒[1, 1]' - 𝐒[1, 2] * 𝐒[1, 2]' - 𝐒[2, 1] * 𝐒[2, 1]' +
           𝐒[2, 2] * 𝐒[2, 2]')
    𝐙₂₃ = -𝐒[1, 1] * 𝐒[1, 2]' + 𝐒[2, 2] * 𝐒[2, 1]'
    𝐙₂₄ = 1.0im * (𝐒[1, 1] * 𝐒[1, 2]' + 𝐒[2, 2] * 𝐒[2, 1]')
    𝐙₃₁ = -𝐒[1, 1] * 𝐒[2, 1]' - 𝐒[2, 2] * 𝐒[1, 2]'
    𝐙₃₂ = -𝐒[1, 1] * 𝐒[2, 1]' + 𝐒[2, 2] * 𝐒[1, 2]'
    𝐙₃₃ = 𝐒[1, 1] * 𝐒[2, 2]' + 𝐒[1, 2] * 𝐒[2, 1]'
    𝐙₃₄ = -1.0im * (𝐒[1, 1] * 𝐒[2, 2]' + 𝐒[2, 1] * 𝐒[1, 2]')
    𝐙₄₁ = 1.0im * (𝐒[2, 1] * 𝐒[1, 1]' + 𝐒[2, 2] * 𝐒[1, 2]')
    𝐙₄₂ = 1.0im * (𝐒[2, 1] * 𝐒[1, 1]' - 𝐒[2, 2] * 𝐒[1, 2]')
    𝐙₄₃ = -1.0im * (𝐒[2, 2] * 𝐒[1, 1]' - 𝐒[1, 2] * 𝐒[2, 1]')
    𝐙₄₄ = 𝐒[2, 2] * 𝐒[1, 1]' - 𝐒[1, 2] * 𝐒[2, 1]'

    𝐙 = @SMatrix [𝐙₁₁ 𝐙₁₂ 𝐙₁₃ 𝐙₁₄; 𝐙₂₁ 𝐙₂₂ 𝐙₂₃ 𝐙₂₄; 𝐙₃₁ 𝐙₃₂ 𝐙₃₃ 𝐙₃₄; 𝐙₄₁ 𝐙₄₂ 𝐙₄₃ 𝐙₄₄]
    return real.(𝐙)
end

@testitem "Can calculate phase matrix from amplitude scattering matrix" begin
    using TransitionMatrices

    @test all(phase_matrix([1+2im 2+3im; 0.2-0.5im 0.5-0.2im]) .≈
              [9.29 -4.0 -8.2 -0.79
               8.71 -4.0 -7.8 -1.21
               0.4 1.2 -1.0 -0.4
               2.8 -1.0 -2.8 1.2])
end

"""
```
albedo(𝐓::AbstractTransitionMatrix)
```

Calculate the single scattering albedo from the given T-Matrix.
"""
function albedo(𝐓::AbstractTransitionMatrix)
    scattering_cross_section(𝐓) / extinction_cross_section(𝐓)
end

@testitem "albedo of" begin
    using TransitionMatrices: Spheroid, transition_matrix, albedo

    @testset "non-absorbing scatterers should be equal to 1.0" begin
        s = Spheroid(1.5, 1.0, complex(1.311))
        𝐓 = transition_matrix(s, 2π, 5, 40)
        @test albedo(𝐓) ≈ 1.0
    end

    @testset "absorbing scatterers should be less than 1.0" begin
        s = Spheroid(1.5, 1.0, 1.5 + 0.01im)
        𝐓 = transition_matrix(s, 2π, 5, 40)
        @test albedo(𝐓) < 1.0
    end
end

"""
```
absorption_cross_section(𝐓::AbstractTransitionMatrix, λ=2π)
```

Calculate the absorption cross section from the given T-Matrix.
"""
function absorption_cross_section(𝐓::AbstractTransitionMatrix, λ = 2π)
    extinction_cross_section(𝐓, λ) - scattering_cross_section(𝐓, λ)
end

@testitem "absorption cross section of" begin
    using TransitionMatrices: Spheroid, transition_matrix, absorption_cross_section

    @testset "non-absorbing scatterers should be approximate to 0.0" begin
        s = Spheroid(1.5, 1.0, complex(1.311))
        𝐓 = transition_matrix(s, 2π, 5, 40)

        # The result may be slightly negative due to numerical errors
        @test abs(absorption_cross_section(𝐓)) < 1e-8
    end

    @testset "absorbing scatterers should be less than 1.0" begin
        s = Spheroid(1.5, 1.0, 1.5 + 0.01im)
        𝐓 = transition_matrix(s, 2π, 5, 40)
        @test absorption_cross_section(𝐓) > 0.0
    end
end

@doc raw"""
```
expansion_coefficients(𝐓::AbstractTransitionMatrix{CT, N}, λ) where {CT, N}
```

Calculate the expansion coefficients from an arbitrary T-Matrix, using Eq. (24) -- (74) in Bi et al. (2014).

Parameters:

- `𝐓`: The precalculated T-Matrix of a scatterer.
- `λ`: The wavelength.

Keyword arguments:

- `full`: Whether to return the full expansion coefficients (`β₃` to `β₆`). Default to `false`.
"""
function expansion_coefficients(𝐓::AbstractTransitionMatrix{CT, N}, λ;
                                full = false) where {CT, N}
    Cˢᶜᵃ = Float64(scattering_cross_section(𝐓, λ))
    λ = Float64(λ)
    ci = OffsetArray([(1im)^(i & 4) for i in (-N):N], (-N):N)
    s = OffsetArray([Float64(2i + 1) for i in 0:(2N)], 0:(2N))
    ss = sqrt.(s)
    sig = OffsetArray([1 - 2 * (i % 2) for i in 0:(4N)], 0:(4N))

    wig_table_init(4N, 3)
    wig_temp_init(4N)

    it = OrderDegreeIterator(N)
    T₁ = OffsetArray(zeros(ComplexF64, 2N + 1, N, 2N + 1, N), (-N):N, 1:N, (-N):N, 1:N)
    T₂ = OffsetArray(zeros(ComplexF64, 2N + 1, N, 2N + 1, N), (-N):N, 1:N, (-N):N, 1:N)
    T₃ = OffsetArray(zeros(ComplexF64, 2N + 1, N, 2N + 1, N), (-N):N, 1:N, (-N):N, 1:N)
    T₄ = OffsetArray(zeros(ComplexF64, 2N + 1, N, 2N + 1, N), (-N):N, 1:N, (-N):N, 1:N)

    for (n, m) in it
        for (n′, m′) in it
            T₁[m, n, m′, n′] = 𝐓[m, n, m′, n′, 1, 1] + 𝐓[m, n, m′, n′, 1, 2] +
                               𝐓[m, n, m′, n′, 2, 1] + 𝐓[m, n, m′, n′, 2, 2]
            T₂[m, n, m′, n′] = 𝐓[m, n, m′, n′, 1, 1] + 𝐓[m, n, m′, n′, 1, 2] -
                               𝐓[m, n, m′, n′, 2, 1] - 𝐓[m, n, m′, n′, 2, 2]
            T₃[m, n, m′, n′] = 𝐓[-m, n, -m′, n′, 1, 1] - 𝐓[-m, n, -m′, n′, 1, 2] -
                               𝐓[-m, n, -m′, n′, 2, 1] + 𝐓[-m, n, -m′, n′, 2, 2]
            T₄[m, n, m′, n′] = 𝐓[-m, n, -m′, n′, 1, 1] - 𝐓[-m, n, -m′, n′, 1, 2] +
                               𝐓[-m, n, -m′, n′, 2, 1] - 𝐓[-m, n, -m′, n′, 2, 2]
        end
    end

    # A₁ = zeros(ComplexF64, N)
    # A₂ = zeros(ComplexF64, N)
    # A₃ = zeros(ComplexF64, N)
    # A₄ = zeros(ComplexF64, N)

    A₁ = OffsetArray(zeros(ComplexF64, 2N + 1, N, N, 2N + 1), (-N):N, 1:N, 1:N, 0:(2N))
    A₂ = OffsetArray(zeros(ComplexF64, 2N + 1, N, N, 2N + 1), (-N):N, 1:N, 1:N, 0:(2N))
    A₃ = OffsetArray(zeros(ComplexF64, 2N + 1, N, N, 2N + 1), (-N):N, 1:N, 1:N, 0:(2N))
    A₄ = OffsetArray(zeros(ComplexF64, 2N + 1, N, N, 2N + 1), (-N):N, 1:N, 1:N, 0:(2N))

    a = [inv(√(2n + 1)) for n in 1:(2N)]

    # #TODO: not finished
    # for n in 1:N
    #     for n₁ in 0:(N + n)
    #         for n′ in max(1, abs(n - n₁)):min(N, n₁ + n)
    #             A₁[n′] = 0.0
    #             A₂[n′] = 0.0
    #             A₃[n′] = 0.0
    #             A₄[n′] = 0.0
    #         end
    #     end
    # end

    for n₁ in 0:(2N), n′ in 1:N, n in 1:N, k in (-N):N
        lb = max(-N - k, -N)
        ub = min(N - k, N)

        A₁[k, n, n′, n₁] = im^((n′ - n) & 3) * a[n′] *
                           sum(clebschgordan(n, m₁, n₁, k, n′) * T₁[m₁, n, m₁ + k, n′]
                               for m₁ in lb:ub)
        A₂[k, n, n′, n₁] = im^((n′ - n) & 3) * a[n′] *
                           sum(clebschgordan(n, m₁, n₁, k, n′) * T₂[m₁, n, m₁ + k, n′]
                               for m₁ in lb:ub)
        A₃[k, n, n′, n₁] = im^((n′ - n) & 3) * a[n′] *
                           sum(clebschgordan(n, m₁, n₁, k, n′) * T₃[m₁, n, m₁ + k, n′]
                               for m₁ in lb:ub)
        A₄[k, n, n′, n₁] = im^((n′ - n) & 3) * a[n′] *
                           sum(clebschgordan(n, m₁, n₁, k, n′) * T₄[m₁, n, m₁ + k, n′]
                               for m₁ in lb:ub)
    end

    B₁ = OffsetArray(zeros(ComplexF64, 2N + 1, 2N + 1, N, 2N + 1), (-N):N, (-N):N, 1:N,
                     0:(2N))
    B₂ = OffsetArray(zeros(ComplexF64, 2N + 1, 2N + 1, N, 2N + 1), (-N):N, (-N):N, 1:N,
                     0:(2N))
    B₃ = OffsetArray(zeros(ComplexF64, 2N + 1, 2N + 1, N, 2N + 1), (-N):N, (-N):N, 1:N,
                     0:(2N))
    B₄ = OffsetArray(zeros(ComplexF64, 2N + 1, 2N + 1, N, 2N + 1), (-N):N, (-N):N, 1:N,
                     0:(2N))

    for n₁ in 0:(2N), n in 1:N, m in (-N):N, k in (-N):N
        lb = max(1, abs(n - n₁))
        ub = min(n + n₁, N)

        B₁[k, m, n, n₁] = sum(clebschgordan(n, m, n₁, 1 - m, n′) * A₁[k, n, n′, n₁]
                              for n′ in lb:ub; init = 0)
        B₂[k, m, n, n₁] = sum(clebschgordan(n, m, n₁, 1 - m, n′) * A₂[k, n, n′, n₁]
                              for n′ in lb:ub; init = 0)
        B₃[k, m, n, n₁] = sum(clebschgordan(n, m, n₁, 1 - m, n′) * A₃[k, n, n′, n₁]
                              for n′ in lb:ub; init = 0)
        B₄[k, m, n, n₁] = sum(clebschgordan(n, m, n₁, 1 - m, n′) * A₄[k, n, n′, n₁]
                              for n′ in lb:ub; init = 0)
    end

    D₀₀ = OffsetArray(zeros(ComplexF64, 2N + 1, N, N), (-N):N, 1:N, 1:N)
    D₀₋₀ = OffsetArray(zeros(ComplexF64, 2N + 1, N, N), (-N):N, 1:N, 1:N)
    D₋₀₋₀ = OffsetArray(zeros(ComplexF64, 2N + 1, N, N), (-N):N, 1:N, 1:N)
    D₂₂ = OffsetArray(zeros(ComplexF64, 2N + 1, N, N), (-N):N, 1:N, 1:N)
    D₂₋₂ = OffsetArray(zeros(ComplexF64, 2N + 1, N, N), (-N):N, 1:N, 1:N)
    D₋₂₋₂ = OffsetArray(zeros(ComplexF64, 2N + 1, N, N), (-N):N, 1:N, 1:N)
    D₀₂ = OffsetArray(zeros(ComplexF64, 2N + 1, N, N), (-N):N, 1:N, 1:N)
    D₋₀₂ = OffsetArray(zeros(ComplexF64, 2N + 1, N, N), (-N):N, 1:N, 1:N)

    for n′ in 1:N, n in 1:N
        for m in (-min(n, n′)):min(n, n′)
            lb = abs(m - 1)
            ub = min(n, n′) + N

            D₀₀[m, n, n′] = sum((2n₁ + 1) *
                                sum(B₃[k, m, n, n₁] * B₃[k, m, n′, n₁]'
                                    for k in max(-N, -n₁):min(N, n₁))
                                for n₁ in lb:ub; init = 0)
            D₀₋₀[m, n, n′] = sum((2n₁ + 1) *
                                 sum(B₂[k, m, n, n₁] * B₂[k, m, n′, n₁]'
                                     for k in max(-N, -n₁):min(N, n₁))
                                 for n₁ in lb:ub; init = 0)
            D₋₀₋₀[m, n, n′] = sum((2n₁ + 1) *
                                  sum(B₁[k, m, n, n₁] * B₁[k, m, n′, n₁]'
                                      for k in max(-N, -n₁):min(N, n₁))
                                  for n₁ in lb:ub; init = 0)
        end

        for m in max(-n, -n′ + 2):min(n, n′ + 2)
            lb = abs(m - 1)
            ub = min(n, n′) + N

            D₂₂[m, n, n′] = sum((2n₁ + 1) *
                                sum(B₁[k, m, n, n₁] * B₃[-k, 2 - m, n′, n₁]'
                                    for k in max(-N, -n₁):min(N, n₁))
                                for n₁ in lb:ub; init = 0)
            D₂₋₂[m, n, n′] = sum((2n₁ + 1) *
                                 sum(B₄[k, m, n, n₁] * B₂[-k, 2 - m, n′, n₁]'
                                     for k in max(-N, -n₁):min(N, n₁))
                                 for n₁ in lb:ub; init = 0)
            D₋₂₋₂[m, n, n′] = sum((2n₁ + 1) *
                                  sum(B₃[k, m, n, n₁] * B₁[-k, 2 - m, n′, n₁]'
                                      for k in max(-N, -n₁):min(N, n₁))
                                  for n₁ in lb:ub; init = 0)
            D₀₂[m, n, n′] = sum((2n₁ + 1) *
                                sum(B₂[k, m, n, n₁] * B₃[-k, 2 - m, n′, n₁]'
                                    for k in max(-N, -n₁):min(N, n₁))
                                for n₁ in lb:ub; init = 0)
            D₋₀₂[m, n, n′] = sum((2n₁ + 1) *
                                 sum(B₁[k, m, n, n₁] * B₄[-k, 2 - m, n′, n₁]'
                                     for k in max(-N, -n₁):min(N, n₁))
                                 for n₁ in lb:ub; init = 0)
        end
    end

    h_const = λ^2 / (Cˢᶜᵃ * 4 * π)
    h = OffsetArray([s[l] * h_const * ss[n] / ss[n′]
                     for l in 0:(2N), n in 1:N, n′ in 1:N],
                    0:(2N),
                    1:N,
                    1:N)

    # Calculate g
    g₀₀ = OffsetArray(zeros(ComplexF64, 2N + 1), 0:(2N))
    g₀₋₀ = OffsetArray(zeros(ComplexF64, 2N + 1), 0:(2N))
    g₋₀₋₀ = OffsetArray(zeros(ComplexF64, 2N + 1), 0:(2N))
    g₂₂ = OffsetArray(zeros(ComplexF64, 2N + 1), 0:(2N))
    g₂₋₂ = OffsetArray(zeros(ComplexF64, 2N + 1), 0:(2N))
    g₋₂₋₂ = OffsetArray(zeros(ComplexF64, 2N + 1), 0:(2N))
    g₀₂ = OffsetArray(zeros(ComplexF64, 2N + 1), 0:(2N))
    g₋₀₂ = OffsetArray(zeros(ComplexF64, 2N + 1), 0:(2N))

    for l in 0:(2N)
        for n in 1:N
            for n′ in max(1, abs(n - l)):(min(n + l, N))
                cg1 = clebschgordan(n, 1, l, 0, n′)
                sm₀₀ = 0.0
                sm₀₋₀ = 0.0
                sm₋₀₋₀ = 0.0

                for m in (-min(n, n′)):min(n, n′)
                    cg = clebschgordan(n, m, l, 0, n′)
                    sm₀₀ += cg * D₀₀[m, n, n′]
                    sm₀₋₀ += cg * D₀₋₀[m, n, n′]
                    sm₋₀₋₀ += cg * D₋₀₋₀[m, n, n′]
                end

                g₀₀[l] += h[l, n, n′] * cg1 * sm₀₀
                g₀₋₀[l] += h[l, n, n′] * cg1 * sig[n + n′ + l] * sm₀₋₀
                g₋₀₋₀[l] += h[l, n, n′] * cg1 * sm₋₀₋₀

                if l ≥ 2
                    cg2 = clebschgordan(n, -1, l, 2, n′)
                    sm₂₂ = 0.0
                    sm₂₋₂ = 0.0
                    sm₋₂₋₂ = 0.0
                    sm₀₂ = 0.0
                    sm₋₀₂ = 0.0

                    for m in max(-n, -n′ + 2):min(n, n′ + 2)
                        cg = clebschgordan(n, -m, l, 2, n′)
                        sm₂₂ += cg * D₂₂[m, n, n′]
                        sm₂₋₂ += cg * D₂₋₂[m, n, n′]
                        sm₋₂₋₂ += cg * D₋₂₋₂[m, n, n′]
                        sm₀₂ += cg * D₀₂[m, n, n′]
                        sm₋₀₂ += cg * D₋₀₂[m, n, n′]
                    end

                    g₂₂[l] += h[l, n, n′] * cg2 * sm₂₂
                    g₂₋₂[l] += h[l, n, n′] * cg2 * sig[n + n′ + l] * sm₂₋₂
                    g₋₂₋₂[l] += h[l, n, n′] * cg2 * sm₋₂₋₂
                    g₀₂[l] += -h[l, n, n′] * cg1 * sm₀₂
                    g₋₀₂[l] += -h[l, n, n′] * cg1 * sig[n + n′ + l] * sm₋₀₂
                end
            end
        end
    end

    α₁ = @. 0.5real(g₀₀ + 2g₀₋₀ + g₋₀₋₀)
    α₂ = @. real(g₂₂ + g₂₋₂)
    α₃ = @. real(g₂₂ - g₂₋₂)
    α₄ = @. 0.5real(g₀₀ - 2g₀₋₀ + g₋₀₋₀)
    β₁ = @. real(g₀₂ + g₋₀₂)
    β₂ = @. imag(g₀₂ - g₋₀₂)
    β₃ = @. -imag(g₀₂ + g₋₀₂)
    β₄ = @. -imag(g₂₂)
    β₅ = @. 0.5real(g₀₀ - g₋₀₋₀)
    β₆ = @. real(g₀₂ - g₋₀₂)

    wig_temp_free()
    wig_table_free()

    if full
        return α₁, α₂, α₃, α₄, β₁, β₂, β₃, β₄, β₅, β₆
    else
        return α₁, α₂, α₃, α₄, β₁, β₂
    end
end

@doc raw"""
```
asymmetry_parameter(𝐓, λ)
```

Calculate the asymmetry parameter from the given transition matrix, using Eq. (4.92) in Mishchenko et al. (2002):

```math
\langle\cos\Theta\rangle=\frac{\alpha_1^1}{3}
```

"""
function asymmetry_parameter(𝐓::AbstractTransitionMatrix, λ)
    α₁, _ = expansion_coefficients(𝐓, λ)
    return α₁[1] / 3
end

@doc raw"""
```
scattering_matrix(𝐓, λ, θs)
```

Calculate expansion coefficients first and then calculate scatterering matrix elements.

Parameters:

- `𝐓`: The transition matrix.
- `λ`: The wavelength.
- `θs`: The scattering angles to be evaluated in degrees.
"""
function scattering_matrix(𝐓::AbstractTransitionMatrix, λ, θs::AbstractVector)
    α₁, α₂, α₃, α₄, β₁, β₂ = expansion_coefficients(𝐓, λ)
    return scattering_matrix(α₁, α₂, α₃, α₄, β₁, β₂, θs)
end

@doc raw"""
```
scattering_matrix(α₁, α₂, α₃, α₄, β₁, β₂, θs)
```

Calculate the scatterering matrix elements from the given expansion coefficients.

Parameters:

- `α₁`, `α₂`, `α₃`, `α₄`, `β₁`, `β₂`: The precalculated expansion coefficients.
- `θs`: The scattering angles to be evaluated in degrees.
"""
function scattering_matrix(α₁, α₂, α₃, α₄, β₁, β₂, θs::AbstractVector)
    lmax = length(α₁) - 1
    θs = deg2rad.(θs)
    Nθ = length(θs)

    F = zeros(Nθ, 6)
    Threads.@threads for i in eachindex(θs)
        θ = θs[i]
        d₀₀ = wigner_d_recursion(0, 0, lmax, θ)
        d₂₂ = wigner_d_recursion(2, 2, lmax, θ)
        d₂₋₂ = wigner_d_recursion(2, -2, lmax, θ)
        d₀₂ = wigner_d_recursion(0, 2, lmax, θ)

        F₁₁ = sum(α₁[l] * d₀₀[l] for l in 0:lmax)
        F₂₂₊₃₃ = sum((α₂[l] + α₃[l]) * d₂₂[l] for l in 2:lmax)
        F₂₂₋₃₃ = sum((α₂[l] - α₃[l]) * d₂₋₂[l] for l in 2:lmax)
        F₂₂ = (F₂₂₊₃₃ + F₂₂₋₃₃) / 2
        F₃₃ = F₂₂₊₃₃ - F₂₂
        F₄₄ = sum(α₄[l] * d₀₀[l] for l in 0:lmax)
        F₁₂ = -sum(β₁[l] * d₀₂[l] for l in 2:lmax)
        F₃₄ = -sum(β₂[l] * d₀₂[l] for l in 2:lmax)

        F[i, :] .= F₁₁, F₁₂, F₂₂, F₃₃, F₃₄, F₄₄
    end

    return F
end

@doc raw"""
```
scattering_matrix(α₁, α₂, α₃, α₄, β₁, β₂, β₃, β₄, β₅, β₆, θs)
```

Calculate all 10 independent scatterering matrix elements (`F₁₁`, `F₁₂`, `F₁₃`, `F₁₄`, `F₂₂`, `F₂₃`, `F₂₄`, `F₃₃`, `F₃₄`, `F₄₄`) from the given expansion coefficients.

The scattering matrix can be expressed as:

```math
\mathbf{F}=\left[\begin{array}{cccc}
F_{11} & F_{12} & F_{13} & F_{14} \\
F_{12} & F_{22} & F_{23} & F_{24} \\
-F_{13} & -F_{23} & F_{33} & F_{34} \\
F_{14} & F_{24} & -F_{34} & F_{44}
\end{array}\right]
```

Parameters:

- `α₁`, `α₂`, `α₃`, `α₄`, `β₁`, `β₂`, `β₃`, `β₄`, `β₅`, `β₆`: The precalculated expansion coefficients.
- `θs`: The scattering angles to be evaluated in degrees.
"""
function scattering_matrix(α₁, α₂, α₃, α₄, β₁, β₂, β₃, β₄, β₅, β₆, θs::AbstractVector)
    lmax = length(α₁) - 1
    θs = deg2rad.(θs)
    Nθ = length(θs)

    F = zeros(Nθ, 10)
    Threads.@threads for i in eachindex(θs)
        θ = θs[i]
        d₀₀ = wigner_d_recursion(0, 0, lmax, θ)
        d₂₂ = wigner_d_recursion(2, 2, lmax, θ)
        d₂₋₂ = wigner_d_recursion(2, -2, lmax, θ)
        d₀₂ = wigner_d_recursion(0, 2, lmax, θ)
        d₂₀ = wigner_d_recursion(2, 0, lmax, θ)

        F₁₁ = sum(α₁[l] * d₀₀[l] for l in 0:lmax)
        F₂₂₊₃₃ = sum((α₂[l] + α₃[l]) * d₂₂[l] for l in 2:lmax)
        F₂₂₋₃₃ = sum((α₂[l] - α₃[l]) * d₂₋₂[l] for l in 2:lmax)
        F₂₂ = (F₂₂₊₃₃ + F₂₂₋₃₃) / 2
        F₃₃ = F₂₂₊₃₃ - F₂₂
        F₄₄ = sum(α₄[l] * d₀₀[l] for l in 0:lmax)
        F₁₂ = -sum(β₁[l] * d₀₂[l] for l in 2:lmax)
        F₃₄ = -sum(β₂[l] * d₀₂[l] for l in 2:lmax)

        F₁₃ = sum(β₃[l] * d₀₂[l] for l in 2:lmax)
        F₂₃ = sum(β₄[l] * d₂₂[l] for l in 2:lmax)
        F₁₄ = sum(β₅[l] * d₀₀[l] for l in 0:lmax)
        F₂₄ = sum(β₆[l] * d₂₀[l] for l in 2:lmax)

        F[i, :] .= F₁₁, F₁₂, F₁₃, F₁₄, F₂₂, F₂₃, F₂₄, F₃₃, F₃₄, F₄₄
    end

    return F
end
