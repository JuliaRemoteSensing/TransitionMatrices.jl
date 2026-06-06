# Examples

These pages are [Pluto.jl](https://plutojl.org) notebooks rendered to HTML (with
plots). Each one is a self-contained walkthrough of one capability; the live,
editable sources live in the [`examples/`](https://github.com/JuliaRemoteSensing/TransitionMatrices.jl/tree/main/examples)
directory of the repository.

## The notebooks

- [**Shapes quick-start gallery**](shapes_gallery.md) — define a spheroid,
  cylinder, Chebyshev particle, and N-fold prism; compute each T-matrix and its
  `Qsca`/`Qext`/`g`/`ω`, with notes on which solver to use.
- [**Solver landscape & convergence**](solver_landscape.md) — the two-layer API
  (`EBCM`, `IITM`, `Iterative`, `stable`) on a high-aspect spheroid, showing
  classic EBCM diverging with `nₘₐₓ` while the stabilized path holds.
- [**Angular scattering & polarization**](angular_scattering.md) — the
  orientation-averaged scattering (Mueller) matrix vs angle (phase function
  `F₁₁`, linear polarization `-F₁₂/F₁₁`, …), plus fixed-geometry amplitude and
  phase matrices.
- [**Orientation averaging**](orientation_averaging.md) — Mishchenko's analytic
  `RandomOrientationTransitionMatrix` vs the numerical `orientation_average`, and
  the numerical average converging to it.
- [**Spectral sensitivity**](spectral_sensitivity.md) — fast wavelength /
  refractive-index sensitivities via the Sh-matrix moment-separation backend:
  `prepare_sh` runs the geometry quadrature once, then `ForwardDiff` over the
  cheap reconstruction gives exact `∂Cₛ𝚌ₐ/∂λ` and `∂Cₛ𝚌ₐ/∂mᵣ` across a spectrum.
- [**Rain radar observables**](rain_radar.md) — a self-contained
  radiative-transfer prototype: spheroidal IITM T-matrices, a rain-drop
  axis-ratio model, water/ice refractive-index fits, and a drop-size distribution
  feeding dual-pol radar moments and brightness-temperature integrals.

## Running them locally

The notebooks use the `examples/` environment, which develops the parent package:

```julia-repl
julia> import Pkg; Pkg.activate("examples"); Pkg.instantiate()  # first time only
julia> using Pluto; Pluto.run()
```

Then open a notebook from the Pluto start page. Each notebook's first cell
activates this environment and develops `TransitionMatrices` from the repository
root, so it works from a fresh checkout.
