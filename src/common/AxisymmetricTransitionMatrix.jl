struct AxisymmetricTransitionMatrix{CT, N, V <: AbstractVector{<:AbstractMatrix{CT}}, T} <:
       AbstractTransitionMatrix{CT, N}
    𝐓::V
end

Base.@propagate_inbounds function Base.getindex(axi::AxisymmetricTransitionMatrix{CT, N, V},
                                                m::Integer, n::Integer, m′::Integer,
                                                n′::Integer, p::Integer,
                                                p′::Integer) where {CT, N, V}
    if m != m′ || abs(m) > min(n, n′)
        zero(CT)
    else
        mₐ = abs(m)
        sig = m >= 0 ? 1 : (-1)^(p + p′)
        nn = N - max(1, mₐ) + 1
        n₁ = (p - 1) * nn + n - max(1, mₐ) + 1
        n₂ = (p′ - 1) * nn + n′ - max(1, mₐ) + 1
        axi.𝐓[mₐ + 1][n₁, n₂] * sig
    end
end
