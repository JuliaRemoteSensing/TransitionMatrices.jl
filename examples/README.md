# Examples

Runnable examples for `TransitionMatrices.jl`, as [Pluto.jl](https://plutojl.org)
notebooks. They use this `examples/` environment ([`Project.toml`](Project.toml)),
which includes `Pluto`, `Plots`, `ForwardDiff`, and develops the parent package.

## Running a notebook

```julia
julia> import Pkg; Pkg.activate("examples"); Pkg.instantiate()  # first time only
julia> using Pluto; Pluto.run()
```

Then open the notebook from the Pluto start page. Each notebook's first cell
activates this environment and develops `TransitionMatrices` from the repo root,
so it works from a fresh checkout (the first run installs dependencies).

## Notebooks

- [`spectral_sensitivity.jl`](spectral_sensitivity.jl) — **fast wavelength /
  refractive-index sensitivities** via the Sh-matrix moment-separation backend.
  `prepare_sh` runs the geometry quadrature once; differentiating the cheap
  reconstruction with `ForwardDiff` then gives exact `∂Cₛ𝚌ₐ/∂λ` and `∂Cₛ𝚌ₐ/∂mᵣ`
  across a whole spectrum at a large speedup over differentiating the
  from-scratch assembly. Includes plots of the spectra and their sensitivities,
  and a timing comparison.
