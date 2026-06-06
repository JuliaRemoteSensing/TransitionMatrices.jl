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

### Prism

The following defines a hexagonal prism (`N = 6` base edges) with base edge
length `a = 1.0`, height `h = 2.0`, and refractive index `m = 1.5 + 0.01im`. The
number of base edges is a type parameter, supplied as the first argument:

```julia
prism = Prism(6, 1.0, 2.0, 1.5 + 0.01im);
```

A prism is **not** axisymmetric, so its T-matrix is built with the [`IITM`](@ref)
solver rather than EBCM (see [Choosing a solver](@ref)).

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
𝐓 = transition_matrix(prism, 2π, IITM(10, 20, 24, 24))   # non-axisymmetric: add Nφ
𝐓 = transition_matrix(spheroid, 2π, ShMatrix(10, 200))   # Sh-matrix moment separation
```

## Choosing a solver

All engines produce the same kind of transition matrix; pick by the scatterer
and the regime:

| Situation | Recommended |
| --- | --- |
| Sphere or coated sphere | Mie (`bhmie` / `bhcoat`) — analytic and exact |
| Axisymmetric, moderate size and aspect ratio | bare `transition_matrix(s, λ)` (auto-converged classic EBCM) |
| High-aspect-ratio **spheroid** | [`Iterative(EBCM; stable = true)`](@ref Iterative) (see the next section) |
| Non-axisymmetric shape (e.g. a prism) | [`IITM`](@ref) — the only built-in route |
| Axisymmetric but EBCM ill-conditioned | [`IITM`](@ref) |
| Many wavelengths / refractive indices for **one** shape | [`ShMatrix`](@ref) / [`prepare_sh`](@ref) — prepare once, reuse |

Two orthogonal knobs control accuracy in hard cases:

- **`stable`** rebuilds the EBCM `𝐔`-matrix integrals with a cancellation-free
  formulation. It is the right tool for high-aspect-ratio spheroids in `Float64`,
  and is **spheroid-only** (see the next section).
- **Extended precision** (`Double64`, `Float128`, `Arb`) raises the size
  parameter EBCM can reach and lowers round-off floors in general. It works for
  any shape but is several times slower. The two stack: `stable` with a
  `Double64` element type reaches `~1e-25` at high aspect ratio.

When in doubt, start with the bare call; if it fails to converge or loses
accuracy, switch to `Iterative(EBCM; stable = true)` for spheroids, or to
[`IITM`](@ref) otherwise.

## High-aspect-ratio spheroids: the `stable` formulation

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

## Parameter sweeps with the Sh-matrix method (`prepare_sh`)

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

After getting the T-Matrix, you can calculate the far-field scattering
properties. Each function has a short `calc_*` alias:

- [`amplitude_matrix`](@ref), a.k.a. `calc_S`
- [`phase_matrix`](@ref), a.k.a. `calc_Z`
- [`extinction_cross_section`](@ref), a.k.a. `calc_Cext`
- [`scattering_cross_section`](@ref), a.k.a. `calc_Csca`
- [`absorption_cross_section`](@ref), a.k.a. `calc_Cabs`
- [`albedo`](@ref), a.k.a. `calc_ω`
- [`asymmetry_parameter`](@ref), a.k.a. `calc_g`
- [`scattering_matrix`](@ref), a.k.a. `calc_F` (orientation-averaged)

The orientation-averaged scalar quantities take the T-matrix and the wavelength:

```julia
𝐓 = transition_matrix(spheroid, 2π)

Csca = scattering_cross_section(𝐓, 2π)
Cext = extinction_cross_section(𝐓, 2π)
Cabs = absorption_cross_section(𝐓, 2π)
ω    = albedo(𝐓)                     # = Csca / Cext
g    = asymmetry_parameter(𝐓, 2π)    # ⟨cos Θ⟩
```

For a **fixed scattering geometry**, the ``2\times2`` amplitude (Jones) matrix is
evaluated at an incidence direction `(ϑᵢ, φᵢ)` and a scattering direction
`(ϑₛ, φₛ)` (all angles in radians); the ``4\times4`` Mueller / phase matrix
follows from it:

```julia
𝐒 = amplitude_matrix(𝐓, 0.0, 0.0, π / 2, 0.0; λ = 2π)   # forward-to-side geometry
𝐙 = phase_matrix(𝐒)
```

The orientation-averaged scattering matrix is evaluated on a grid of scattering
angles ``\Theta``:

```julia
θs = range(0, π; length = 181)
𝐅 = scattering_matrix(𝐓, 2π, θs)     # one matrix per angle
```

## Near-field reconstruction

Beyond the far field, the **electromagnetic field around the particle** can be
reconstructed from any T-matrix. For an incident plane wave the scattered field
is expanded in radiating vector spherical wave functions with coefficients
``(p, q) = \mathbf{T}\,(a, b)``, where ``(a, b)`` are the plane-wave coefficients;
the total external field is ``\mathbf{E}_\text{inc} + \mathbf{E}_\text{sca}``.

The incident plane wave propagates along ``\hat{\mathbf n} = (\vartheta_\text{inc},
\varphi_\text{inc})`` with a Jones polarization `(Eθ, Eφ)` in the spherical basis
``(\hat{\boldsymbol\vartheta}, \hat{\boldsymbol\varphi})`` at that direction. Field
points are Cartesian (any `AbstractVector` of length 3); the polar axis is `+z`.

```julia
𝐓 = transition_matrix(spheroid, 2π)
λ = 2π

# +z propagation, x-polarized (θ̂ at ϑ=0 is x̂):
ϑ_inc, φ_inc, Eθ, Eφ = 0.0, 0.0, 1.0 + 0im, 0.0im
r⃗ = [4.0, 0.0, 3.0]                                   # a point outside the particle

E_inc = incident_field(λ, ϑ_inc, φ_inc, Eθ, Eφ, r⃗)   # analytic plane wave
E_sca = scattered_field(𝐓, λ, ϑ_inc, φ_inc, Eθ, Eφ, r⃗)
E_tot = total_field(𝐓, λ, ϑ_inc, φ_inc, Eθ, Eφ, r⃗)   # = E_inc + E_sca
```

When evaluating a **dense grid** of field points, compute the scattered
coefficients once and reuse them:

```julia
p, q = scattering_coefficients(𝐓, ϑ_inc, φ_inc, Eθ, Eφ)
field = [scattered_field(p, q, λ, [x, 0.0, z]) for z in zs, x in xs]
```

The underlying VSWFs are available directly as [`vswf`](@ref) / [`vswf_cartesian`](@ref).

!!! warning "Region of validity"
    The radiating (scattered) expansion converges only **outside the smallest
    sphere circumscribing the particle** (the Rayleigh hypothesis). For a sphere
    that boundary is the surface; for a non-spherical particle, evaluate only
    outside its circumscribing sphere.

For the field **inside** a homogeneous sphere, [`internal_field`](@ref)
reconstructs it from the analytic Mie internal coefficients (at the internal
wavenumber ``k_\text{int} = m_r k``); it is continuous, tangentially, with the
external [`total_field`](@ref) across the surface:

```julia
x, mᵣ = 3.0, 1.5 + 0.05im                        # size parameter and index
r⃗ = [0.3, 0.2, 0.5]                              # a point inside the sphere

E_in = internal_field(x, mᵣ, λ, ϑ_inc, φ_inc, Eθ, Eφ, r⃗)

# reuse the coefficients across a dense interior grid:
c, d = internal_coefficients(x, mᵣ, ϑ_inc, φ_inc, Eθ, Eφ)
E_in = internal_field(c, d, mᵣ, λ, r⃗)
```

For a general **axisymmetric** particle the internal field is reconstructed from
the EBCM matrices (``\mathbf{c} = \tfrac{1}{2}\mathbf{Q}^{-1}\mathbf{a}`` per
azimuthal block), via the shape-based methods:

```julia
spheroid = Spheroid(1.0, 2.0, 1.4 + 0.02im)
E_in = internal_field(spheroid, λ, nmax, Ng, ϑ_inc, φ_inc, Eθ, Eφ, r⃗)
# or, reused across a grid:
c, d = internal_coefficients(spheroid, λ, nmax, Ng, ϑ_inc, φ_inc, Eθ, Eφ)
E_in = internal_field(c, d, spheroid.m, λ, r⃗)
```

!!! note "Region of validity and conditioning"
    The interior expansion is mathematically guaranteed within the inscribed sphere
    (radius `rmin(shape)`); empirically it stays stable and physical throughout the
    interior (it does **not** diverge past the inscribed sphere). The binding limit
    is instead EBCM ``\mathbf{Q}``-matrix conditioning at high aspect ratio with
    high `nmax` — since ``\mathbf{c} = \tfrac12\mathbf{Q}^{-1}\mathbf{a}`` inherits
    ``\mathbf{Q}^{-1}``, the same Float64 cancellation that limits the T-matrix
    corrupts the coefficients everywhere — so keep `nmax`/aspect where the T-matrix
    is reliable. The convention is validated against Mie on the degenerate sphere
    (`Spheroid(R, R, mᵣ)`); absolute accuracy beyond the inscribed sphere for
    non-spherical shapes is not independently verified.

See the [Near-field maps from a T-matrix](examples/near_field.md) example for a
field-enhancement map.

## Orientation averaging

For randomly oriented particles you usually want orientation-averaged
quantities. There are two routes:

```julia
𝐓 = transition_matrix(spheroid, 2π)

# analytic random-orientation average — fast and exact (Mishchenko et al. (2002), Eq. 5.96)
𝐓ᵣ = RandomOrientationTransitionMatrix(𝐓)

# numerical average over an orientation distribution p(α, β, γ); uniform here
𝐓ₙ = orientation_average(𝐓, (α, β, γ) -> 1 / (8π^2); Nα = 40, Nβ = 40, Nγ = 1)
```

The numerical average converges to the analytic one as the angular grid is
refined (see the *Orientation averaging* example in [Examples](@ref)). The cross
sections and `scattering_matrix` above apply to either result.

A single **fixed** orientation is obtained by rotating the T-matrix by Euler
angles (Z-Y-Z), without recomputing it from the shape (`RotZYZ` is re-exported):

```julia
𝐓rot = rotate(𝐓, RotZYZ(0.3, 0.5, 0.0))
𝐒    = amplitude_matrix(𝐓rot, 0.0, 0.0, π / 2, 0.0)
```

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
