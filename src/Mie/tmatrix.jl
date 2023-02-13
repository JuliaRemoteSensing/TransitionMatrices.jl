struct MieTransitionMatrix{CT, N, V <: AbstractVector{CT}} <: AbstractTransitionMatrix{CT, N}
    a::V
    b::V
end

function MieTransitionMatrix{CT, N}(x, m) where {CT, N}
    T = real(CT)
    a, b = bhmie(T, x, m; nₘₐₓ = N)
    MieTransitionMatrix{CT, N, Vector{CT}}(a, b)
end 

Base.copy(𝐓::MieTransitionMatrix{CT, N}) where {CT, N} = MieTransitionMatrix{CT, N, Vector{CT}}(copy(𝐓.a), copy(𝐓.b))
Base.@propagate_inbounds function Base.getindex(mie::MieTransitionMatrix{CT, N, V}, m::Integer, n::Integer, m′::Integer, n′::Integer, p::Integer, p′::Integer) where {CT, N, V}
    if m != m′ || n != n′ || p != p′ || abs(m) > n
        zero(CT)
    else
        p == 1 ? -mie.b[n] : -mie.a[n]
    end
end

rotate(𝐓::MieTransitionMatrix, ::Rotation{3}) = copy(𝐓)

@testitem "MieTransitionMatrix" begin
    using TransitionMatrices: MieTransitionMatrix, RotZYZ, TransitionMatrix, rotate

    @testset "remains the same under rotations" begin
        𝐓 = MieTransitionMatrix{ComplexF64, 5}(1.0, 1.311)
        @test all(isapprox.(𝐓, rotate(TransitionMatrix{ComplexF64, 5, typeof(𝐓)}(𝐓), RotZYZ(0.2, 0.3, 0.4)); atol=eps(Float64)))
        @test all(isapprox.(𝐓, rotate(TransitionMatrix{ComplexF64, 5, typeof(𝐓)}(𝐓), RotZYZ(0.8, 0.0, -1.0)); atol=eps(Float64)))
    end
end
