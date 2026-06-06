### A Pluto.jl notebook ###
# v1.0.1

using Markdown
using InteractiveUtils

# ╔═╡ 2173870e-6082-11f1-89cd-4976913f2274
begin
    import Pkg
    Pkg.activate(@__DIR__)
    Pkg.develop(; path = dirname(@__DIR__))
    Pkg.instantiate()
    using TransitionMatrices, Plots
end

# ╔═╡ 21738632-6082-11f1-8c1a-07c117ac9474
md"""
# Shapes quick-start gallery

`TransitionMatrices.jl` handles several scatterer families. This gallery defines
one of each, computes its T-matrix, and reports the basic far-field efficiencies
``Q_\text{sca}``, ``Q_\text{ext}``, the asymmetry parameter ``g`` and the
single-scattering albedo ``\omega``.

The axisymmetric shapes (spheroid, cylinder, Chebyshev particle) use the default
EBCM solver via `calc_T`; the N-fold **prism** is not axisymmetric, so it uses
the IITM solver.
"""

# ╔═╡ 21738718-6082-11f1-b6a5-5d3fd6d47663
md"""
## Define the shapes

All at wavelength `λ = 2π` (so the size parameter equals the radius) and the same
refractive index `m = 1.5 + 0.02im`.
"""

# ╔═╡ 2173872c-6082-11f1-808e-0fe676b65624
begin
    λ = 2π
    m = 1.5 + 0.02im
    shapes = (
        spheroid = TransitionMatrices.Spheroid{Float64, ComplexF64}(2.0, 1.0, m),
        cylinder = TransitionMatrices.Cylinder{Float64, ComplexF64}(1.0, 2.0, m),
        chebyshev = TransitionMatrices.Chebyshev{Float64, ComplexF64}(1.0, 0.1, 4, m),
        prism6 = TransitionMatrices.Prism(6, 1.0, 2.0, m)
    )
end

# ╔═╡ 21738736-6082-11f1-a54e-010814dd76d9
md"""
## Compute and tabulate

`calc_T` auto-converges the axisymmetric shapes; the prism uses a fixed-order
IITM solve. Pluto renders the resulting vector of named tuples as a table.
"""

# ╔═╡ 2173874a-6082-11f1-8f4b-f3cee77b5fae
begin
    function Tmatrix(s)
        s isa TransitionMatrices.Prism ?
        transition_matrix(s, λ, IITM(8, 12, 24, 24)) : calc_T(s, λ)
    end

    summary = map(collect(pairs(shapes))) do (name, s)
        T = Tmatrix(s)
        (; shape = name,
            Qsca = round(calc_Csca(T, λ); digits = 4),
            Qext = round(calc_Cext(T, λ); digits = 4),
            g = round(asymmetry_parameter(T, λ); digits = 4),
            ω = round(albedo(T); digits = 4))
    end
end

# ╔═╡ 21738752-6082-11f1-8f31-fb7e99c75cc0
md"""
## Compare the efficiencies
"""

# ╔═╡ 2173875e-6082-11f1-8933-7da031c63642
let
    names = String.(getindex.(summary, :shape))
    bar(names, [getindex.(summary, :Qsca) getindex.(summary, :Qext)];
        label = ["Qsca" "Qext"], ylabel = "efficiency", legend = :topright,
        title = "far-field efficiencies by shape", size = (680, 380))
end

# ╔═╡ 21738768-6082-11f1-8536-4125de072774
md"""
## Which solver?

- **Spheroid / Cylinder / Chebyshev** — axisymmetric: EBCM (`calc_T`). For
  high-aspect spheroids add the stabilized path, `Iterative(EBCM; stable=true)`.
- **Prism / arbitrary 3-D** — not axisymmetric: IITM (`IITM(nₘₐₓ, Nr, Nϑ, Nφ)`).
- For wavelength / refractive-index **sweeps** of an axisymmetric shape, see the
  Sh-matrix examples (`ShMatrix`, `prepare_sh`).
"""

# ╔═╡ Cell order:
# ╟─21738632-6082-11f1-8c1a-07c117ac9474
# ╠═2173870e-6082-11f1-89cd-4976913f2274
# ╟─21738718-6082-11f1-b6a5-5d3fd6d47663
# ╠═2173872c-6082-11f1-808e-0fe676b65624
# ╟─21738736-6082-11f1-a54e-010814dd76d9
# ╠═2173874a-6082-11f1-8f4b-f3cee77b5fae
# ╟─21738752-6082-11f1-8f31-fb7e99c75cc0
# ╠═2173875e-6082-11f1-8933-7da031c63642
# ╟─21738768-6082-11f1-8536-4125de072774
