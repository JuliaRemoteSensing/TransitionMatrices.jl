"""
Iterator for the order-degree pairs of the given maximum order `nₘₐₓ`.
"""
struct OrderDegreeIterator
    nₘₐₓ::Int
end

Base.iterate(::OrderDegreeIterator) = ((1, -1), (1, -1))
function Base.iterate(iter::OrderDegreeIterator, (n, m))
    if m == n
        n == iter.nₘₐₓ ? nothing : ((n + 1, -n - 1), (n + 1, -n - 1))
    else
        ((n, m + 1), (n, m + 1))
    end
end

Base.firstindex(iter::OrderDegreeIterator) = 1
Base.lastindex(iter::OrderDegreeIterator) = length(iter)
function Base.getindex(iter::OrderDegreeIterator, idx)
    n = floor(Int, √idx)
    m = idx - n^2 - n
    (n, m)
end
Base.length(iter::OrderDegreeIterator) = iter.nₘₐₓ * (iter.nₘₐₓ + 2)
Base.size(x::OrderDegreeIterator) = (length(x),)
Base.eltype(::OrderDegreeIterator) = Tuple{Int, Int}
Base.isdone(iter::OrderDegreeIterator, state) = state >= (iter.nₘₐₓ, iter.nₘₐₓ)

@testitem "OrderDegreeIterator" begin
    using TransitionMatrices: OrderDegreeIterator

    @test iterate(OrderDegreeIterator(3)) == ((1, -1), (1, -1))
    @test iterate(OrderDegreeIterator(3), (1, 1)) == ((2, -2), (2, -2))
    @test collect(OrderDegreeIterator(2)) ==
          [(1, -1), (1, 0), (1, 1), (2, -2), (2, -1), (2, 0), (2, 1), (2, 2)]
    @test size(OrderDegreeIterator(100)) == (10200,)
    @test !Base.isdone(OrderDegreeIterator(2), (1, 1))
    @test Base.isdone(OrderDegreeIterator(2), (2, 2))
end

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

- For a `MieTransitionMatrix`, the underlying Mie coefficients are copied and a new `MieTransitionMatrix` will be returned.
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
                    sign = iseven(m′ + m₂) ? 1 : -1
                    𝐓′[m, n, m′, n′, p, p′] += coeff[m, m₁] * d[m, m₁, n] *
                                               conj(coeff[m′, m₂]) * d[m₂, m′, n′] * sign *
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
