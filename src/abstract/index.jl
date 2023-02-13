@doc raw"""
A general T-Matrix ``T_{m n m^{\prime} n^{\prime}}^{k l}`` stored in a 6-dimensional array, in the order ``(m, n, m^{\prime}, n^{\prime}, k, l)``.
"""
abstract type AbstractTransitionMatrix{T, N} <: AbstractArray{T, 6} end

order(::AbstractTransitionMatrix{T, N}) where {T, N} = N

@doc raw"""
Rotate a general T-Matrix using Eq. (5.29) of Mishchenko et al. (2002).

```math
T_{m n m^{\prime} n^{\prime}}^{k l}(L ; \alpha, \beta, \gamma)=\sum_{m_1=-n}^n \sum_{m_2=-n^{\prime}}^{n^{\prime}} D_{m m_1}^n(\alpha, \beta, \gamma) T_{m_1 n m_2 n^{\prime}}^{k l}(P) D_{m_2 m^{\prime}}^{n^{\prime}}(-\gamma,-\beta,-\alpha)\quad k,l=1,2
```
"""
function rotate(𝐓::AbstractTransitionMatrix{C, N}, rot::Rotation{3, F}) where {C, F, N}
    # Get the Euler angle in Z-Y-Z order.
    zyz = RotZYZ(rot)
    α, β, γ = zyz.theta1, zyz.theta2, zyz.theta3

    # Calculate the wigner-d functions that will be used.
    d = OffsetArray(zeros(C, 2N + 1, 2N + 1, N + 1), -N:N, -N:N, 0:N)
    for m in -N:N
        for m′ in -N:N
            sₘᵢₙ = max(abs(m), abs(m′))
            wigner_d_recursion!(view(d, m, m′, sₘᵢₙ:N), m, m′, N, β)
        end
    end

    # Calculate the coefficients used for wigner-D functions
    coeff = OffsetArray([cis(-(m * α + m′ * γ)) for m in -N:N, m′ in -N:N], -N:N, -N:N)

    # Calculate the rotated T-Matrix
    𝐓′ = similar(𝐓)
    fill!(𝐓′, 0)
    for k in 1:2, l in 1:2
        for n′ in 0:N
            for m′ in -n′:n′
                for n in 0:N
                    for m in -n:n
                        tmp = 0.0
                        for m₂ in -n′:n′, m₁ in -n:n
                            tmp += coeff[m, m₁] * d[m, m₁, n] * 𝐓′[m₁, n, m₂, n′, k, l] * conj(coeff[m′, m₂]) * d[m₂, m′, n′]
                        end
                        𝐓′[m, n, m′, n′, k, l] = tmp
                    end
                end
            end
        end
    end

    𝐓′
end
