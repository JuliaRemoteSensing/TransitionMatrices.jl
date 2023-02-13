@doc raw"""
Calculate Wigner (small) d-function ``d_{mn}^s(\theta)`` for a single ``(m, n, s)`` combination, using Eq. (B.1) of Mishchenko et al. (2002).

```math
\begin{aligned}
d_{m n}^{s}(\vartheta)=& \sqrt{(s+m) !(s-m) !(s+n) !(s-n) !} \\
& \times \sum_{k=\max(0,m-n)}^{\min(s + m, s - n)}(-1)^{k} \frac{\left(\cos \frac{1}{2} \vartheta\right)^{2 s-2 k+m-n}\left(\sin \frac{1}{2} \vartheta\right)^{2 k-m+n}}{k !(s+m-k) !(s-n-k) !(n-m+k) !}
\end{aligned}
```
"""
function wigner_d(::Type{T}, m::Integer, n::Integer, s::Integer, θ::Number) where {T}
    if s < max(abs(m), abs(n))
        return zero(T)
    end

    θ = T(θ)
    half = θ / 2
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

@inline function wigner_d(m::Integer, n::Integer, s::Integer, θ::Number)
    wigner_d(Float64, m, n, s, θ)
end

@doc raw"""
Calculate Wigner (small) d-function ``d_{mn}^s(\theta)`` for ``s\in[s_{\min}=\max(|m|, |n|),s_{\max}]`` (and also its derivative) via upward recursion, using Eq. (B.22) of Mishchenko et al. (2002).

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
function wigner_d_recursion(::Type{T}, m::Integer, n::Integer, smax::Integer, θ::Number;
                            deriv::Bool = false) where {T}
    smax >= max(abs(m), abs(n)) || error("Error: smax < max(|m|, |n|)")

    θ = T(θ)
    cosθ = cos(θ)
    sinθ = sin(θ)
    smin = max(abs(m), abs(n))
    d = zeros(T, smax - smin + 1)
    if deriv
        d_deriv = zeros(T, smax - smin + 1)
    end

    d = OffsetArray(d, smin:smax)
    if deriv
        d_deriv = OffsetArray(d_deriv, smin:smax)
    end

    ξ = m <= n ? 1 : (-1)^((m - n) & 1)
    d₀ = 0
    d₁ = d[smin] = ξ * T(2)^(-smin) *
                   √(factorial(T, 2smin) / factorial(T, abs(m - n)) /
                     factorial(T, abs(m + n))) *
                   (1 - cosθ)^(abs(m - n) / 2) * (1 + cosθ)^(abs(m + n) / 2)

    for s in smin:smax
        s1m = √T((s + 1)^2 - m^2)
        s1n = √T((s + 1)^2 - n^2)
        sm = √T(s^2 - m^2)
        sn = √T(s^2 - n^2)

        if s == 0
            d₂ = T(2s + 1) * cosθ * d₁
        else
            d₂ = 1 / (s * s1m * s1n) *
                 ((2s + 1) * (s * (s + 1) * cosθ - m * n) * d₁ -
                  (s + 1) * sm * sn * d₀)
        end
        if s < smax
            d[s + 1] = d₂
        end

        if deriv
            if s == 0
                d_deriv[s] = 0
            else
                d_deriv[s] = 1 / sinθ * (-(s + 1) * sm * sn / (s * (2s + 1)) * d₀ -
                              T(m * n // (s * (s + 1))) * d₁ +
                              s * s1m * s1n / ((s + 1) * (2s + 1)) * d₂)
            end
        end

        d₀, d₁ = d₁, d₂
    end

    return deriv ? (d, d_deriv) : d
end

@inline function wigner_d_recursion(m::Integer, n::Integer, smax::Integer, θ::Number;
                                    deriv::Bool = false)
    return wigner_d_recursion(Float64, m, n, smax, θ; deriv)
end

function wigner_d_recursion!(d::AbstractVector{T}, m::Integer, n::Integer, smax::Integer,
                             θ::Number) where {T}
    smax >= max(abs(m), abs(n)) || error("Error: smax < max(|m|, |n|)")

    θ = T(θ)
    cosθ = cos(θ)
    smin = max(abs(m), abs(n))
    d = OffsetArray(d, smin:smax)

    ξ = m <= n ? 1 : (-1)^((m - n) & 1)
    d₀ = 0
    d₁ = d[smin] = ξ * T(2)^(-smin) *
                   √(factorial(T, 2smin) / factorial(T, abs(m - n)) /
                     factorial(T, abs(m + n))) *
                   (1 - cosθ)^(abs(m - n) / 2) * (1 + cosθ)^(abs(m + n) / 2)

    for s in smin:smax
        s1m = √T((s + 1)^2 - m^2)
        s1n = √T((s + 1)^2 - n^2)
        sm = √T(s^2 - m^2)
        sn = √T(s^2 - n^2)

        if s == 0
            d₂ = T(2s + 1) * cosθ * d₁
        else
            d₂ = 1 / (s * s1m * s1n) *
                 ((2s + 1) * (s * (s + 1) * cosθ - m * n) * d₁ -
                  (s + 1) * sm * sn * d₀)
        end
        if s < smax
            d[s + 1] = d₂
        end

        d₀, d₁ = d₁, d₂
    end

    return d
end

@doc raw"""
Wigner D-function ``D^j_{mn}(\theta)`` defined as:

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
    α = T(α)
    β = T(β)
    γ = T(γ)
    return cis(-(m * α + m′ * γ)) * wigner_d(T, m, m′, n, β)
end

@inline function wigner_D(m::Integer, m′::Integer, n::Integer, α::Number, β::Number,
                          γ::Number)
    return wigner_D(Float64, m, m′, n, α, β, γ)
end

function wigner_D_recursion(::Type{T}, m::Integer, m′::Integer, nmax::Integer, α::Number,
                            β::Number,
                            γ::Number) where {T}
    α = T(α)
    β = T(β)
    γ = T(γ)
    d = wigner_d_recursion(T, m, m′, nmax, β)
    return cis(-(m * α + m′ * γ)) * d
end

@inline function wigner_D_recursion(m::Integer, m′::Integer, nmax::Integer, α::Number,
                                    β::Number, γ::Number)
    return wigner_D_recursion(Float64, m, m′, nmax, α, β, γ)
end

function wigner_D_recursion!(d::AbstractVector{T}, m::Integer, m′::Integer, nmax::Integer,
                             α::Number, β::Number,
                             γ::Number) where {T}
    α = T(α)
    β = T(β)
    γ = T(γ)
    wigner_d_recursion!(d, m, m′, nmax, β)
    factor = cis(-(m * α + m′ * γ))
    d .*= factor
    return d
end
