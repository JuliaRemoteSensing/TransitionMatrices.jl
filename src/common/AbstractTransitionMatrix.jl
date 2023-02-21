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
Rotate the given T-Matrix `𝐓` by the Euler angle `rot` and generate a new T-Matrix.

### General T-Matrix

```
rotate(𝐓::AbstractTransitionMatrix{CT, N}, rot::Rotation{3})
```

For a general T-Matrix, Eq. (5.29) in Mishchenko et al. (2002) is used as a fallback. A `TransitionMatrix` will be returned, which is the most general yet concrete type.

```math
T_{m n m^{\prime} n^{\prime}}^{p p′}(L ; \alpha, \beta, \gamma)=\sum_{m_1=-n}^n \sum_{m_2=-n^{\prime}}^{n^{\prime}} D_{m m_1}^n(\alpha, \beta, \gamma) T_{m_1 n m_2 n^{\prime}}^{p p′}(P) D_{m_2 m^{\prime}}^{n^{\prime}}(-\gamma,-\beta,-\alpha)\quad p,p′=1,2
```

### Axisymmetric T-Matrix

### Mie T-Matrix

```
rotate(𝐓::MieTransitionMatrix{CT, N}, rot::Rotation{3})
```

For a `MieTransitionMatrix`, the underlying Mie coefficients are copied and a new `MieTransitionMatrix` will be returned.

### Orientation-averaged T-Matrix

```
rotate(oa::OrientationAveragedTransitionMatrix{CT, N, V}, ::Rotation{3}) where {CT, N, V} =
    typeof(oa)(copy(oa.𝐓))
```

For an `OrientationAveragedTransitionMatrix`, the underlying T-Matrix is copied and a new `OrientationAveragedTransitionMatrix` will be returned.

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
Calculate the amplitude matrix of the given T-Matrix `𝐓` at the given incidence and scattering angles. `k₁` is the wavenumber of the incident wave in the host medium, which should be calculated by `k₁ = 2π * mₕ / λ`, where `mₕ` is the refractive index of the host medium and `λ` is the wavelength of the incident wave. The default value is `k₁ = 1.0`.

### General T-Matrix

```
amplitude_matrix(𝐓::AbstractTransitionMatrix{CT, N}, ϑᵢ, φᵢ, ϑₛ, φₛ, k₁=1.0)
```

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

### Axisymmetric T-Matrix

### Mie T-Matrix

"""
function amplitude_matrix(𝐓::AbstractTransitionMatrix{CT, N}, ϑᵢ, φᵢ, ϑₛ, φₛ,
                          k₁ = 1.0) where {CT, N}
    T = real(CT)
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
Calculate the orientation average of a transition matrix, given the orientation distribution function ``p_o(\alpha,\beta,\gamma)``. 

### General T-Matrix

```
orientation_average(𝐓::AbstractTransitionMatrix{CT, N}, pₒ; Nα = 10, Nβ = 10, Nγ = 10) where {CT, N}
```

For a general T-Matrix and a general orientation distribution function, we use numerical integration to calculate the orientation average. The orientation average is given by

```math
\langle T_{m n m^{\prime} n^{\prime}}^{p p^{\prime}}(L)\rangle = \int_0^{2\pi}\mathrm{d}\alpha\int_0^{\pi}\mathrm{d}\beta\sin\beta\int_0^{2\pi}\mathrm{d}\gamma p_o(\alpha,\beta,\gamma) T_{m n m^{\prime} n^{\prime}}^{p p^{\prime}}(L; \alpha,\beta,\gamma)
```

Parameters:

- `𝐓`: The T-Matrix to be orientation averaged.
- `pₒ`: The orientation distribution function. Note that the ``\sin\beta`` part is already included.
- `Nα`: The number of points used in the numerical integration of ``\alpha``. Default to 10.
- `Nβ`: The number of points used in the numerical integration of ``\beta``. Default to 10.
- `Nγ`: The number of points used in the numerical integration of ``\gamma``. Default to 10.

!!! note

    This is the fallback method and does not utilize any symmetry, so it is expected to be slow. You should use specified versions of this function, or implement your own if there is no suited version for your combination of T-Matrix and orientation distribution function.

    You may also need to test the convergence of `Nα`, `Nβ` and `Nγ` manually. If any one is too small, there will be large errors in the results.

### Random orientation T-Matrix and Mie T-Matrix

```
orientation_average(mie::MieTransitionMatrix, _pₒ; _kwargs...)
orientation_average(mie::MieTransitionMatrix, _pₒ; _kwargs...)
```

Both types are invariant under rotation. Therefore, the original T-Matrix will be returned.
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
