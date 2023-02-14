@doc raw"""
`wigner_d([T=Float64,], m::Integer, n::Integer, s::Integer, ϑ::Number) where {T}`

Calculate Wigner (small) d-function ``d_{mn}^s(\theta)`` for a single ``(m, n, s)`` combination, using Eq. (B.1) of Mishchenko et al. (2002).

```math
\begin{aligned}
d_{m n}^{s}(\vartheta)=& \sqrt{(s+m) !(s-m) !(s+n) !(s-n) !} \\
& \times \sum_{k=\max(0,m-n)}^{\min(s + m, s - n)}(-1)^{k} \frac{\left(\cos \frac{1}{2} \vartheta\right)^{2 s-2 k+m-n}\left(\sin \frac{1}{2} \vartheta\right)^{2 k-m+n}}{k !(s+m-k) !(s-n-k) !(n-m+k) !}
\end{aligned}
```
"""
function wigner_d(::Type{T}, m::Integer, n::Integer, s::Integer, ϑ::Number) where {T}
    if s < max(abs(m), abs(n))
        return zero(T)
    end

    ϑ = T(ϑ)
    half = ϑ / 2
    cos_half = cos(half)
    sin_half = sin(half)
    kmin = max(0, m - n)
    kmax = min(s + m, s - n)
    sig = (-1)^(kmin & 1)
    ch = cos_half^(2s - 2kmin + m - n)
    sh = sin_half^(2kmin - m + n)
    fac = factorial(T, kmin) * factorial(T, s + m - kmin) * factorial(T, s - n - kmin) *
          factorial(T, n - m + kmin)
    d = sig * ch * sh / fac
    for k in (kmin + 1):kmax
        sig = -sig
        ch *= cos_half^(-2)
        sh *= sin_half^2
        fac *= k * (n - m + k)
        fac /= (s + m - k + 1) * (s - n - k + 1)
        d += sig * ch * sh / fac
    end

    return d * √(factorial(T, s + m) * factorial(T, s - m) * factorial(T, s + n) *
             factorial(T, s - n))
end

@inline function wigner_d(m::Integer, n::Integer, s::Integer, ϑ::Number)
    wigner_d(Float64, m, n, s, ϑ)
end

@doc raw"""
`wigner_d_recursion([T=Float64,], m::Integer, n::Integer, smax::Integer, ϑ::Number; deriv::Bool = false) where {T}`

Calculate Wigner (small) d-function ``d_{mn}^s(\theta)`` for ``s\in[s_{\min}=\max(|m|, |n|),s_{\max}]`` (alternatively, its derivative as well) via upward recursion, using Eq. (B.22) of Mishchenko et al. (2002).

```math
\begin{aligned}
d_{m n}^{s+1}(\vartheta)=& \frac{1}{s \sqrt{(s+1)^{2}-m^{2}} \sqrt{(s+1)^{2}-n^{2}}}\left\{(2 s+1)[s(s+1) x-m n] d_{m n}^{s}(\vartheta)\right.\\
&\left.-(s+1) \sqrt{s^{2}-m^{2}} \sqrt{s^{2}-n^{2}} d_{m n}^{s-1}(\vartheta)\right\}, \quad s \geq s_{\min }
\end{aligned}
```

The initial terms are given by Eq. (B.23) and Eq. (B.24).

```math
\begin{array}{l}
d_{m n}^{s_{\min }-1}(\vartheta)=0 \\
d_{m n}^{s_{\min }}(\vartheta)=\xi_{m n} 2^{-s_{\min }}\left[\frac{\left(2 s_{\min }\right) !}{(|m-n|) !(|m+n|) !}\right]^{1 / 2}(1-x)^{|m-n| / 2}(1+x)^{|m+n| / 2}
\end{array}
```

where

```math
\xi_{m n}=\left\{\begin{array}{ll}
1 & \text { for } n \geq m \\
(-1)^{m-n} & \text { for } n<m
\end{array}\right.
```
"""
function wigner_d_recursion(::Type{T}, m::Integer, n::Integer, smax::Integer, ϑ::Number;
                            deriv::Bool = false) where {T}
    smax >= max(abs(m), abs(n)) || error("Error: smax < max(|m|, |n|)")
    smin = max(abs(m), abs(n))

    d = zeros(T, smax - smin + 1)
    d_deriv = deriv ? zeros(T, smax - smin + 1) : nothing

    wigner_d_recursion!(d, m, n, smax, ϑ; deriv = d_deriv)
end

@inline function wigner_d_recursion(m::Integer, n::Integer, smax::Integer, ϑ::Number;
                                    deriv::Bool = false)
    return wigner_d_recursion(Float64, m, n, smax, ϑ; deriv)
end

"""
`wigner_d_recursion!(d::AbstractVector{T}, m::Integer, n::Integer, smax::Integer, ϑ::Number; deriv=nothing)`

Calculate the Wigner d-function recursively, in place.
"""
function wigner_d_recursion!(d::AbstractVector{T}, m::Integer, n::Integer, smax::Integer,
                             ϑ::Number; deriv = nothing) where {T}
    smax >= max(abs(m), abs(n)) || error("Error: smax < max(|m|, |n|)")

    ϑ = T(ϑ)
    cosϑ = cos(ϑ)
    !isnothing(deriv) && abs(cosϑ) ≈ 1 && min(abs(m), abs(n)) != 0 && error("Error: wigner_d_recursion! can only calculate derivatives for |cosϑ| < 1 when |m| != 0 and |n| != 0")

    sinϑ = sin(ϑ)
    smin = max(abs(m), abs(n))
    d = OffsetArray(d, smin:smax)
    if !isnothing(deriv)
        deriv = OffsetArray(deriv, smin:smax)
    end

    ξ = m <= n ? 1 : (-1)^((m - n) & 1)
    d₀ = 0
    d₁ = d[smin] = ξ * T(2)^(-smin) *
                   √(factorial(T, 2smin) / factorial(T, abs(m - n)) /
                     factorial(T, abs(m + n))) *
                   (1 - cosϑ)^(abs(m - n) / 2) * (1 + cosϑ)^(abs(m + n) / 2)

    for s in smin:smax
        s1m = √T((s + 1)^2 - m^2)
        s1n = √T((s + 1)^2 - n^2)
        sm = √T(s^2 - m^2)
        sn = √T(s^2 - n^2)

        if s == 0
            d₂ = T(2s + 1) * cosϑ * d₁
        else
            d₂ = 1 / (s * s1m * s1n) *
                 ((2s + 1) * (s * (s + 1) * cosϑ - m * n) * d₁ -
                  (s + 1) * sm * sn * d₀)
        end
        if s < smax
            d[s + 1] = d₂
        end

        if !isnothing(deriv)
            if s == 0
                deriv[s] = 0
            elseif abs(cosϑ) ≈ 1
                deriv[s] = abs(m) + abs(n) == 1 ? √T(s * (s + 1)) / 2 : 0
            else
                deriv[s] = 1 / sinϑ * (-(s + 1) * sm * sn / (s * (2s + 1)) * d₀ -
                            T(m * n // (s * (s + 1))) * d₁ +
                            s * s1m * s1n / ((s + 1) * (2s + 1)) * d₂)
            end
        end

        d₀, d₁ = d₁, d₂
    end

    return isnothing(deriv) ? d : (d, deriv)
end

@testitem "Wigner d-function" begin
    using TransitionMatrices: wigner_d, wigner_d_recursion, wigner_d_recursion!

    @testset "D($m, 0, n, 0.2) is correct" for m in 0:1
        d₁ = [wigner_d(m, 0, n, 0.2) for n in max(0, m):5]
        d₂ = wigner_d_recursion(m, 0, 5, 0.2)
        @test all(d₁ .≈ collect(d₂))

        d₃ = zeros(length(d₂))
        wigner_d_recursion!(d₃, m, 0, 5, 0.2)
        @test all(d₁ .≈ d₃)
    end
end

@doc raw"""
Calculate the Wigner D-function ``D^j_{mn}(\theta)``, which is defined as:

```math
D_{m m^{\prime}}^{n}(\alpha, \beta, \gamma)=\mathrm{e}^{-\mathrm{i} m \alpha} d_{m m^{\prime}}^{n}(\beta) \mathrm{e}^{-\mathrm{i} m^{\prime} \gamma}
```

where

```math
0 \leq \alpha<2 \pi, \quad 0 \leq \beta \leq \pi, \quad 0 \leq \gamma<2 \pi
```
"""
@inline function wigner_D(::Type{T}, m::Integer, m′::Integer, n::Integer, α::Number,
                          β::Number, γ::Number) where {T}
    return cis(-(m * T(α) + m′ * T(γ))) * wigner_d(T, m, m′, n, β)
end

@inline function wigner_D(m::Integer, m′::Integer, n::Integer, α::Number, β::Number,
                          γ::Number)
    return wigner_D(Float64, m, m′, n, α, β, γ)
end

"""
`wigner_D_recursion([T=Float64,], m::Integer, m′::Integer, nmax::Integer, α::Number, β::Number, γ::Number)`

Calculate the Wigner D-function recursively (use `wigner_d_recursion`).
"""
function wigner_D_recursion(::Type{T}, m::Integer, m′::Integer, nmax::Integer, α::Number,
                            β::Number,
                            γ::Number) where {T}
    d = wigner_d_recursion(T, m, m′, nmax, β)
    return cis(-(m * T(α) + m′ * T(γ))) * d
end

@inline function wigner_D_recursion(m::Integer, m′::Integer, nmax::Integer, α::Number,
                                    β::Number, γ::Number)
    return wigner_D_recursion(Float64, m, m′, nmax, α, β, γ)
end

"""
`wigner_D_recursion!(d::AbstractVector{CT}, m::Integer, m′::Integer, nmax::Integer, α::Number, β::Number, γ::Number)`

Calculate the Wigner D-function recursively, in place.
"""
function wigner_D_recursion!(d::AbstractVector{CT}, m::Integer, m′::Integer, nmax::Integer,
                             α::Number, β::Number,
                             γ::Number) where {CT}
    T = real(CT)
    α = T(α)
    β = T(β)
    γ = T(γ)
    wigner_d_recursion!(d, m, m′, nmax, β)
    factor = cis(-(m * α + m′ * γ))
    d .*= factor
    return d
end

@testitem "Wigner D-function" begin
    using TransitionMatrices: wigner_D, wigner_D_recursion, wigner_D_recursion!

    @testset "D($m, 0, n, 0.2, 0.3, 0.4) is correct" for m in 0:1
        D₁ = [wigner_D(m, 0, n, 0.2, 0.3, 0.4) for n in max(0, m):5]
        D₂ = wigner_D_recursion(m, 0, 5, 0.2, 0.3, 0.4)
        @test all(D₁ .≈ collect(D₂))

        D₃ = zeros(ComplexF64, length(D₂))
        wigner_D_recursion!(D₃, m, 0, 5, 0.2, 0.3, 0.4)
        @test all(D₁ .≈ D₃)
    end
end

@doc raw"""
`pi_func([T=Float64,], m::Integer, n::Integer, ϑ::Number; d=nothing)`

Calculate

```math
\pi_{m n}(\vartheta)=\frac{m d_{0 m}^n(\vartheta)}{\sin \vartheta}
```

If `d` is given, it is used as the value of ``d_{0 m}^n(\vartheta)``.
"""
function pi_func(::Type{T}, m::Integer, n::Integer, ϑ::Number; d = nothing) where {T}
    ϑ = T(ϑ)
    cosϑ = cos(ϑ)
    if cosϑ ≈ 1
        return abs(m) == 1 ? √T(n * (n + 1)) / 2 : zero(T)
    elseif cosϑ ≈ -1
        return abs(m) == 1 ? -√T(n * (n + 1)) / 2 : zero(T)
    else
        if isnothing(d)
            d = wigner_d_recursion(T, 0, m, n, ϑ)[n]
        end
        return m * d / sin(ϑ)
    end
end

@inline function pi_func(m::Integer, n::Integer, ϑ::Number; d = nothing)
    return pi_func(Float64, m, n, ϑ; d = d)
end

@doc raw"""
`tau_func([T=Float64,], m::Integer, n::Integer, ϑ::Number)`

Calculate

```math
\tau_{m n}(\vartheta)=\frac{\mathrm{d} d_{0 m}^n(\vartheta)}{\mathrm{d} \vartheta}
```
"""
function tau_func(::Type{T}, m::Integer, n::Integer, ϑ::Number) where {T}
    _, d′ = wigner_d_recursion(T, 0, m, n, ϑ; deriv=true)
    return d′[n]
end

@inline function tau_func(m::Integer, n::Integer, ϑ::Number)
    return tau_func(Float64, m, n, ϑ)
end
