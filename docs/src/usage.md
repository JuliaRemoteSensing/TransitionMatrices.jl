```@meta
CurrentModule = TransitionMatrices
```

# Usage

## Define a shape

### Spheroid

The following defines a spheroid with semi-major axis `a=1.0`, semi-minor axis `c=0.5`, and a refractive index of `m=1.311`.

```julia
spheroid = Spheroid{Float64, ComplexF64}(1.0, 0.5, 1.311);
```

The type specification `{Float64, ComplexF64}` can be omitted, and in that case, the real and complex types used in the calculations will be inferred from the input parameters.

!!! note

    In the example above, if `{Float64, ComplexF64}` is ommited, do remember to use `complex(1.311)` instead of 1.311 to denote that the refractive index is a complex value. Otherwise errors would occur in later calculations.

### Cylinder

The following defines a cylinder with radius `r=1.0`, height `h=2.0`, and a refractive index of `m=1.5+0.01im`. Note that the type specification `{Float128, ComplexF128}` is used to specify the real and complex types used in the calculations.

```julia
cylinder = Cylinder{Float128, ComplexF128}(1.0, 2.0, 1.5+0.01im);
```

### Chebyshev particle

The following defines a Chebyshev particle with radius `r₀=1.0`, deformation coefficient `ε=0.1`, order `n=4`, and a refractive index of `m=1.5+0.01im`. Note that the type specification `{Arb, Acb}` is used to specify the real and complex types used in the calculations.

```julia
chebyshev = Chebyshev{Arb, Acb}(1.0, 0.1, 4, 1.5+0.01im);
```

## Calculate the T-Matrix

After defining a shape, the T-Matrix is computed with `transition_matrix` (alias
`calc_T`). The solvers are organized into **two layers**, selected by the kind of
*solver object* you pass as the third argument:

- an **iterative solver** ([`Iterative`](@ref)) sweeps the discretization until
  the scattering/extinction efficiencies converge;
- a **fixed solution** ([`EBCM`](@ref), [`IITM`](@ref), [`ShMatrix`](@ref))
  builds the T-matrix once at an explicit discretization.

The bare call uses an automatically-converged classic EBCM — the common default:

```julia
𝐓 = calc_T(spheroid, 2π)                          # ≡ transition_matrix(spheroid, 2π, Iterative(EBCM))
```

To see the convergence process, set `ENV["JULIA_DEBUG"] = "TransitionMatrices"`
before running. Tune the convergence (or switch on stabilization) by passing an
explicit `Iterative`:

```julia
𝐓 = transition_matrix(spheroid, 2π, Iterative(EBCM; threshold = 1e-7))
```

For a fixed build, pass a fixed solver with its discretization:

```julia
𝐓 = transition_matrix(spheroid, 2π, EBCM(10, 100))       # order 10, 100 quadrature points
𝐓 = transition_matrix(cylinder, 2π, IITM(10, 20, 24))    # IITM: order, Nr, Nϑ
𝐓 = transition_matrix(spheroid, 2π, ShMatrix(10, 200))   # Sh-matrix moment separation
```

### High-aspect-ratio spheroids (`stable`)

For spheroids of high aspect ratio the standard EBCM surface integrals lose all
precision in `Float64`: the irregular Riccati–Bessel products develop large
Laurent terms that cancel analytically on integration but not numerically. The
`stable` option assembles the `𝐔`-matrix with the cancellation-free `F⁺`
formulation of [Somerville, Auguié & Le Ru (2013)](https://doi.org/10.1016/j.jqsrt.2012.07.017),
which removes the cancellation and restores a relative accuracy of about `1e-9`
in `Float64`, independent of aspect ratio. It is available on **both layers** —
auto-converged or fixed:

```julia
prolate = Spheroid{Float64, ComplexF64}(2.5198421, 10.079368, 1.55 + 0.01im)  # aspect 4
𝐓 = transition_matrix(prolate, 2π, Iterative(EBCM; stable = true))  # converged + stabilized
𝐓 = transition_matrix(prolate, 2π, EBCM(32, 400; stable = true))    # fixed + stabilized
```

`stable` is **only valid for `Spheroid`** (the cancellation relies on the spheroid
surface) and costs roughly 2–3× the default assembly, so it is opt-in. For a fixed
build to help, choose `nₘₐₓ` large enough that `nₘₐₓ + 1 ≳ k·c + 15`, where `c` is
the largest semi-axis — comparable to the order needed for convergence at that
size anyway (the iterative form finds it for you).

`stable` and extended precision are **orthogonal** and stack. Just raising
precision without `stable` (e.g. `Spheroid{Double64}`) only buys ~15 extra digits
of headroom against the cancellation — enough for moderate aspect ratios but not
high ones — and is several times slower; for a high-aspect spheroid `stable` in
`Float64` is both more accurate and faster. When you need more than `~1e-9`,
combine the two: `stable` with a `Double64` element type reaches `~1e-25` at high
aspect ratio, and also lowers the residual `Float64` round-off floor that appears
as the refractive index approaches `s → 1`.

### Parameter sweeps with the Sh-matrix method (`prepare_sh`)

When you need the T-matrix of one fixed shape at *many* wavelengths or
refractive indices — a spectrum, a dispersion curve, a refractive-index scan —
the bulk of the work (the Gauss–Legendre surface quadrature and the Wigner
functions) does not change between points. The Sh-matrix moment-separation
method (Farafonov family) exploits this: every EBCM surface integral factors as

```math
\int = \sum_{\text{terms}} \big[\text{coefficient in } (k, m_r)\big]\times\big[\text{shape-only moment } M(q)\big],
```

where the *shape moments* depend only on the particle geometry and the azimuthal
index — not on the wavelength or the refractive index. `prepare_sh` computes the
moments once; each `transition_matrix(prep, λ, mᵣ)` is then a cheap
coefficient×moment reconstruction:

```julia
spheroid = Spheroid{Float64, ComplexF64}(2.0, 1.0, 1.5 + 0.02im)
prep = prepare_sh(spheroid, 10, 200)

# one wavelength / refractive index
𝐓 = transition_matrix(prep, 2π, 1.5 + 0.02im)

# a full sweep, reusing the single preparation
λs = range(2π, 3π; length = 50)
𝐓s = transition_matrix_spectrum(prep, λs, 1.5 + 0.02im)        # fixed index
𝐓s = transition_matrix_spectrum(prep, λs, my_dispersion.(λs))  # λ-dependent index
```

The same sweep is available through the [`ShMatrix`](@ref) fixed-solver façade,
which prepares once internally and reconstructs at each point:

```julia
𝐓s = transition_matrix(spheroid, λs, my_dispersion.(λs), ShMatrix(10, 200))
```

For a spheroid the reconstruction also inherits the high-aspect-ratio
stabilization for free: the negative-`r`-power moments vanish analytically, so —
accumulated at high precision inside `prepare_sh` — they contribute ≈0 and the
irregular-product cancellation that defeats the standard `Float64` assembly
never forms (`prep.stable` reports this). The cross sections of the
reconstructed T-matrix match the standard assembly to roughly `1e-9`–`1e-15` for
moderate sizes, and the high-aspect path tracks the `stable=true` result.

`prepare_sh` accepts a keyword `B` (the number of radial power-series terms,
default `max(30, nₘₐₓ+15)`). Like the `F⁺` evaluation, the series needs more
terms as the size parameter `k·rₘₐₓ·|mᵣ|` grows, so the method is best suited to
small-to-moderate size parameters; very large particles still want the recursion
path. The analytic `𝐔` stabilization is enabled only for `Spheroid` — for other
axisymmetric shapes the moment machinery and `𝐏` are correct and the sweep reuse
still applies, but the `𝐔` reconstruction is not stabilized.

## Post-processing

After getting the T-Matrix, you can calculate the far-field scattering properties using the following functions:

- [`amplitude_matrix`](@ref), a.k.a. `calc_S`
- [`phase_matrix`](@ref), a.k.a. `calc_Z`
- [`extinction_cross_section`](@ref), a.k.a. `calc_Cext`
- [`scattering_cross_section`](@ref), a.k.a. `calc_Csca`
- [`absorption_cross_section`](@ref), a.k.a. `calc_Cabs`
- [`albedo`](@ref), a.k.a. `calc_ω`
- [`asymmetry_parameter`](@ref), a.k.a. `calc_g`

And the orientation-averaged scattering matrix:

- [`scattering_matrix`](@ref), a.k.a. `calc_F`

## Differentiation

`TransitionMatrices.jl` supports numerical automatic differentiation using
[`ForwardDiff.jl`](https://github.com/JuliaDiff/ForwardDiff.jl), and also
provides analytical linearization backends for supported solver and parameter
combinations.

The analytical path is exposed through `LinearizationProblem`,
`linearize_transition_matrix`, and `linearize_observable`. It currently covers
Mie size/material/wavelength variables, fixed-order EBCM slices for spheroids,
cylinders, and Chebyshev particles, and fixed-geometry IITM material/wavelength
slices. Unsupported combinations throw `UnsupportedLinearization` instead of
silently falling back to finite differences. For backend support boundaries, see
[Linearization Framework](@ref).

For example, Mie transition-matrix derivatives can be requested directly:

```julia
problem = LinearizationProblem(
    [1.7, 1.311, 0.02, 2π];
    variables = (:x, :mᵣ, :mᵢ, :λ),
) do x
    (; x = x[1], m = complex(x[2], x[3]), λ = x[4], nₘₐₓ = 5)
end

result = linearize_transition_matrix(problem, MieLinearization())
∂T_∂λ = derivative(result, :λ)
```

Here is an example of calculating the gradient of the scattering cross section with respect to `a`, `c`, `mᵣ`, `mᵢ` and `λ`:

```jldoctest
using ForwardDiff, TransitionMatrices

function f(x)
    s = Spheroid(x[1], x[2], complex(x[3], x[4]))
    T₀ = TransitionMatrices.transition_matrix(s, x[5], 10, 100)
    Csca = calc_Csca(T₀)
end

gradient = round.(ForwardDiff.gradient(f, [2.0, 3.0, 1.311, 0.02, 2π]); digits = 12)

# output

5-element Vector{Float64}:
  19.872009666583
   5.740718904092
  89.130506909932
 -45.845687597164
  -9.066448506675
```
