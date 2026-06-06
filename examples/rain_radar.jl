### A Pluto.jl notebook ###
# v1.0.1

using Markdown
using InteractiveUtils

# ╔═╡ 14fa0b92-232b-4a39-b509-88df2453872a
begin
    import Pkg
    Pkg.activate(@__DIR__)
    Pkg.develop(; path = dirname(@__DIR__))
    Pkg.instantiate()
    using TransitionMatrices, Plots
end

# ╔═╡ 65279bf2-9dd8-4a9f-837d-2b15670cc301
md"""
# Rain radar observables

This notebook turns single-particle T-matrices into simple rain radar and
brightness-temperature observables. It is adapted from the radiative-transfer
prototype contributed by
[`@xiongyuup`](https://github.com/xiongyuup) in
[`PR #3`](https://github.com/JuliaRemoteSensing/TransitionMatrices.jl/pull/3).

The original contribution demonstrates a larger workflow with WRF input,
hydrometeor microphysics, radar moments, and microwave radiative transfer. This
example keeps only the self-contained core: water refractive index, rain-drop
axis ratio, spheroidal T-matrices, a Marshall-Palmer-type drop-size
distribution, and integrated dual-polarization moments.
"""

# ╔═╡ 96b07bb0-d7f7-4fdb-9c58-7e43a6a0bf31
md"""
## Water, ice, and drop models

The dielectric fits below are compact microwave models for liquid water and
ice. Rain drops are represented as oblate spheroids with diameter ``D`` and
axis ratio ``c/a``.
"""

# ╔═╡ 109c0c9d-967c-489c-846e-9194a0b6e8a2
begin
    c0_m_per_s() = 2.99792458e8
    c0_mm_per_s() = 1000 * c0_m_per_s()

    function water_refractive_index(freq_hz, T_celsius)
        λ_m = c0_m_per_s() / freq_hz
        ε0 = 8.854e-12
        εs = 78.54 * (1 - 4.579e-3 * (T_celsius - 25) +
              1.19e-5 * (T_celsius - 25)^2 -
              2.8e-8 * (T_celsius - 25)^3)
        ε∞ = 5.27137 + 2.16474e-2 * T_celsius - 1.31198e-3 * T_celsius^2
        α = -16.8129 / (T_celsius + 273) + 6.09265e-2
        λs = 3.3836e-6 * exp(2513.98 / (T_celsius + 273))
        σ = 1.1117e-4

        x = (λs / λ_m)^(1 - α)
        denominator = 1 + 2x * sinpi(α / 2) + (λs / λ_m)^(2 - 2α)
        εreal = ε∞ + (εs - ε∞) * (1 + x * sinpi(α / 2)) / denominator
        εimag = (εs - ε∞) * x * cospi(α / 2) / denominator +
                σ * λ_m / (2π * c0_m_per_s() * ε0)

        sqrt(complex(εreal, εimag))
    end

    function ice_refractive_index(freq_hz, T_celsius)
        λ_m = c0_m_per_s() / freq_hz
        ε0 = 8.854e-12
        εs = 203.168 + 2.5 * T_celsius + 0.15T_celsius^2
        ε∞ = 3.168
        α = 0.288 + 5.2e-3 * T_celsius + 2.3e-4 * T_celsius^2
        λs = 9.990288e-6 * exp(6643.5 / (T_celsius + 273))
        σ = 1.1146e-13 * exp(-6291.2 / (T_celsius + 273))

        x = (λs / λ_m)^(1 - α)
        denominator = 1 + 2x * sinpi(α / 2) + (λs / λ_m)^(2 - 2α)
        εreal = ε∞ + (εs - ε∞) * (1 + x * sinpi(α / 2)) / denominator
        εimag = (εs - ε∞) * x * cospi(α / 2) / denominator +
                σ * λ_m / (2π * c0_m_per_s() * ε0)

        sqrt(complex(εreal, εimag))
    end

    function rain_axis_ratio(D_mm)
        D_mm < 0.7 && return 1.0
        D_mm < 1.5 && return 1.173 - 0.5265D_mm + 0.4698D_mm^2 -
                             0.1317D_mm^3 - 8.5e-3D_mm^4
        1.065 - 6.25e-2D_mm - 3.99e-3D_mm^2 +
        7.66e-4D_mm^3 - 4.095e-5D_mm^4
    end
end

# ╔═╡ 07e1ca0a-3bf8-45a9-ae4d-8c51971c6f57
begin
    round_complex(z; digits = 4) = complex(round(real(z); digits), round(imag(z); digits))

    material_summary = [
        (; material = "water", freq_GHz = 9.375, T = 0.0,
            m = round_complex(water_refractive_index(9.375e9, 0.0))),
        (; material = "ice", freq_GHz = 9.375, T = -10.0,
            m = round_complex(ice_refractive_index(9.375e9, -10.0)))
    ]
end

# ╔═╡ 5d8211a5-3b7a-49a3-8249-474d7a5dc833
md"""
## Single-particle scattering table

Each diameter is solved as a spheroid with the IITM backend. The backward
amplitudes feed reflectivity and differential reflectivity; the forward
amplitudes feed differential phase; extinction feeds the brightness-temperature
toy calculation.
"""

# ╔═╡ 0c68e61c-e0d7-43ea-9c06-f144d042aa9b
begin
    function scattering_table(Ds, freq_hz; temperature_celsius = 0.0,
            solver = IITM(6, 10, 16))
        λ = c0_mm_per_s() / freq_hz
        m = water_refractive_index(freq_hz, temperature_celsius)
        map(Ds) do D
            a = D / 2
            axis_ratio = rain_axis_ratio(D)
            c = a * axis_ratio
            shape = Spheroid{Float64, ComplexF64}(a, c, ComplexF64(m))
            T = transition_matrix(shape, λ, solver)
            Sback = amplitude_matrix(T, π / 2, 0.0, π / 2, π; λ)
            Sfwd = amplitude_matrix(T, π / 2, 0.0, π / 2, 0.0; λ)
            (; D_mm = D,
                axis_ratio,
                Shh_back = Sback[2, 2],
                Svv_back = Sback[1, 1],
                Shh_forward = Sfwd[2, 2],
                Svv_forward = Sfwd[1, 1],
                Cext_mm2 = calc_Cext(T, λ))
        end
    end

    spacing(table) = length(table) > 1 ? table[2].D_mm - table[1].D_mm : 1.0
end

# ╔═╡ 8fce0f52-2fac-4798-84f8-75d084e27d85
begin
    freq_hz = 9.375e9
    temperature_celsius = 0.0
    Ds = collect(0.75:0.75:3.75)
    solver = IITM(6, 10, 16)

    table = scattering_table(Ds, freq_hz; temperature_celsius, solver)
    table_summary = [(; D = row.D_mm,
                         axis_ratio = round(row.axis_ratio; digits = 3),
                         Cext = round(row.Cext_mm2; digits = 4))
                     for row in table]
end

# ╔═╡ aaef0a83-7657-4c40-a875-bc7f18a7b42e
table_summary

# ╔═╡ d56cbf67-bad5-45aa-a2b4-2aa4e184fa8d
md"""
## Bulk rain observables

For compactness this uses
``N(D, R) = N_0 \exp[-\Lambda(R)D]`` with ``D`` in millimeters. The constants are
chosen to produce a plausible monotonic example, not a complete retrieval
algorithm.
"""

# ╔═╡ d78d1d83-508b-4521-ae44-9068d15f4d14
begin
    function drop_distribution(D_mm, rain_rate_mm_h)
        rain_rate_mm_h == 0 ? 0.0 : 8.0e3 * exp(-4.1 * rain_rate_mm_h^(-0.21) * D_mm)
    end

    function dualpol_moments(table, freq_hz, rain_rate_mm_h;
            temperature_celsius = 0.0)
        λ = c0_mm_per_s() / freq_hz
        ΔD = spacing(table)
        m = water_refractive_index(freq_hz, temperature_celsius)
        Kw2 = abs2((m^2 - 1) / (m^2 + 2))
        weights = [drop_distribution(row.D_mm, rain_rate_mm_h) * ΔD
                   for row in table]

        hh = sum(abs2(row.Shh_back) * weights[i] for (i, row) in pairs(table))
        vv = sum(abs2(row.Svv_back) * weights[i] for (i, row) in pairs(table))
        Zhh = hh * λ^4 / (π^5 * Kw2)
        Zvv = vv * λ^4 / (π^5 * Kw2)
        cov = sum(conj(row.Shh_back) * row.Svv_back * weights[i]
        for (i, row) in pairs(table))
        ρhv = abs(cov) / sqrt(hh * vv)
        kdp = 180 / π * 1e-3 * λ *
              sum(real(row.Shh_forward - row.Svv_forward) * weights[i]
              for (i, row) in pairs(table))

        (; R = rain_rate_mm_h,
            ZH = round(10log10(Zhh); digits = 2),
            ZDR = round(10log10(Zhh / Zvv); digits = 2),
            ρhv = round(ρhv; digits = 4),
            KDP = round(kdp; digits = 4))
    end

    function brightness_temperature(table, rain_rate_mm_h;
            path_km = 5.0, Tatm_K = 270.0)
        ΔD = spacing(table)
        α = sum(drop_distribution(row.D_mm, rain_rate_mm_h) *
                row.Cext_mm2 * ΔD for row in table) / 1e6
        τ = α * path_km
        Tatm_K * (1 - exp(-τ))
    end
end

# ╔═╡ 383e1078-e216-4810-b83d-cd7a6590749b
begin
    rain_rates = [1.0, 10.0, 50.0]
    radar = [dualpol_moments(table, freq_hz, R; temperature_celsius)
             for R in rain_rates]
    tb = [(; R, TB = round(brightness_temperature(table, R); digits = 4))
          for R in rain_rates]
    (; radar, tb)
end

# ╔═╡ 38b8df2b-2f47-4f81-b027-f46f5e6e4b75
md"""
## Multi-frequency brightness-temperature lookup

The larger PR #3 prototype writes lookup tables to MAT files. Here the same
idea stays in memory: build one scattering table per frequency, integrate over
the rain distribution, and return a compact table that Pluto renders directly.
"""

# ╔═╡ d6d595e6-4e3f-4c4d-9412-032aef0d940e
begin
    function brightness_temperature_lut(Ds, frequencies_hz, rain_rates;
            temperature_celsius = 0.0, solver = IITM(6, 10, 16),
            path_km = 5.0, Tatm_K = 270.0)
        rows = NamedTuple[]
        for f in frequencies_hz
            local_table = scattering_table(Ds, f; temperature_celsius, solver)
            for R in rain_rates
                push!(rows,
                    (; freq_GHz = round(f / 1e9; digits = 3), R,
                        TB = round(
                            brightness_temperature(local_table, R;
                                path_km, Tatm_K); digits = 4)))
            end
        end
        rows
    end

    lookup_frequencies = [10e9, 18e9]
    tb_lut = brightness_temperature_lut(Ds, lookup_frequencies, rain_rates;
        temperature_celsius, solver)
end

# ╔═╡ 20ae5289-72e6-4f64-b4e2-b61b6eb62473
tb_lut

# ╔═╡ 88f653ee-6d22-4568-8285-57acf2d6c644
md"""
## Plots
"""

# ╔═╡ a16f870a-5a48-4638-8491-524afc5f265f
let
    p1 = plot(getfield.(table_summary, :D), getfield.(table_summary, :Cext);
        marker = :circle, lw = 2, label = "Cext",
        xlabel = "drop diameter D (mm)", ylabel = "extinction cross section (mm²)")
    p2 = plot(getfield.(radar, :R), getfield.(radar, :ZH);
        marker = :circle, lw = 2, xscale = :log10, label = "ZH",
        xlabel = "rain rate (mm h⁻¹)", ylabel = "reflectivity (dBZ)")
    p3 = plot(getfield.(radar, :R), getfield.(radar, :ZDR);
        marker = :circle, lw = 2, xscale = :log10, label = "ZDR",
        xlabel = "rain rate (mm h⁻¹)", ylabel = "differential reflectivity (dB)")
    p4 = plot(getfield.(tb, :R), getfield.(tb, :TB);
        marker = :circle, lw = 2, xscale = :log10, label = "TB",
        xlabel = "rain rate (mm h⁻¹)", ylabel = "brightness temperature (K)")
    plot(p1, p2, p3, p4; layout = (2, 2), size = (820, 620))
end

# ╔═╡ dcf844c1-6d42-47e3-9760-a4d094646796
let
    p = plot(; xscale = :log10, xlabel = "rain rate (mm h⁻¹)",
        ylabel = "brightness temperature (K)", legend = :topleft)
    for f in getfield.(tb_lut, :freq_GHz) |> unique
        subset = filter(row -> row.freq_GHz == f, tb_lut)
        plot!(p, getfield.(subset, :R), getfield.(subset, :TB);
            marker = :circle, lw = 2, label = "$(f) GHz")
    end
    p
end

# ╔═╡ 52c0c938-0a88-41a9-bab8-ed9abf42adee
md"""
## References

- H. Li, Y. Xiong, and Y. Chen, "Simulation of Complex Meteorological Target
  Echoes for Airborne Dual-Polarization Weather Radar Based on Invariant
  Imbedding T-Matrix," *IEEE Transactions on Geoscience and Remote Sensing*,
  62, 5105817, 2024.
- G. Zhang, *Weather Radar Polarimetry*, China Meteorological Press, 2018,
  pp. 39-40.
- B. R. Brown, M. M. Bell, and A. J. Frambach, "Validation of Simulated
  Hurricane Drop Size Distributions Using Polarimetric Radar," *Geophysical
  Research Letters*, 42, 2016.
"""

# ╔═╡ Cell order:
# ╟─65279bf2-9dd8-4a9f-837d-2b15670cc301
# ╠═14fa0b92-232b-4a39-b509-88df2453872a
# ╟─96b07bb0-d7f7-4fdb-9c58-7e43a6a0bf31
# ╠═109c0c9d-967c-489c-846e-9194a0b6e8a2
# ╠═07e1ca0a-3bf8-45a9-ae4d-8c51971c6f57
# ╟─5d8211a5-3b7a-49a3-8249-474d7a5dc833
# ╠═0c68e61c-e0d7-43ea-9c06-f144d042aa9b
# ╠═8fce0f52-2fac-4798-84f8-75d084e27d85
# ╠═aaef0a83-7657-4c40-a875-bc7f18a7b42e
# ╟─d56cbf67-bad5-45aa-a2b4-2aa4e184fa8d
# ╠═d78d1d83-508b-4521-ae44-9068d15f4d14
# ╠═383e1078-e216-4810-b83d-cd7a6590749b
# ╟─38b8df2b-2f47-4f81-b027-f46f5e6e4b75
# ╠═d6d595e6-4e3f-4c4d-9412-032aef0d940e
# ╠═20ae5289-72e6-4f64-b4e2-b61b6eb62473
# ╟─88f653ee-6d22-4568-8285-57acf2d6c644
# ╠═a16f870a-5a48-4638-8491-524afc5f265f
# ╠═dcf844c1-6d42-47e3-9760-a4d094646796
# ╟─52c0c938-0a88-41a9-bab8-ed9abf42adee
