#!/usr/bin/env julia

using NCDatasets
using Serialization
using SpecialFunctions
using TransitionMatrices

include(joinpath(@__DIR__, "shared.jl"))
using .RadiativeTransferShared: MASS_COEFF, RAIN_N0, R_DRY_AIR, WATER_DENSITY,
                                air_density, write_var

const C0_MM_PER_S = 2.99792458e11
const C0_M_PER_S = 2.99792458e8
const FILL = Float32(-999)

function usage()
    println("""
    Compute dual-pol rain radar fields from a WRF-like NetCDF file.

    Usage:
      julia --project=examples/radiative_transfer \\
        examples/radiative_transfer/wrf_dualpol_radar.jl --input fake_wrfout.nc [options]

    Options:
      --input PATH       Input WRF-like NetCDF path
      --output PATH      Output NetCDF path [dualpol_radar.nc]
      --dsd MODE         DSD mode: single or double [double]
      --freq HZ          Radar frequency in Hz [9.375e9]
      --diameters SPEC   Diameter grid start:step:stop in mm [0.75:0.75:3.75]
      --cache PATH       Serialized scattering-table cache [scattering_table.jls]
      --nmax N           IITM order [6]
      --nr N             IITM radial grid [10]
      --ntheta N         IITM polar grid [16]
      --help             Show this message
    """)
end

function parse_args(args)
    opts = Dict(
        "input" => "",
        "output" => "dualpol_radar.nc",
        "dsd" => "double",
        "freq" => "9.375e9",
        "diameters" => "0.75:0.75:3.75",
        "cache" => "scattering_table.jls",
        "nmax" => "6",
        "nr" => "10",
        "ntheta" => "16")
    i = 1
    while i <= length(args)
        arg = args[i]
        if arg == "--help" || arg == "-h"
            usage()
            exit(0)
        elseif startswith(arg, "--") && haskey(opts, arg[3:end])
            i == length(args) && error("missing value for $arg")
            opts[arg[3:end]] = args[i + 1]
            i += 2
        else
            error("unknown argument: $arg")
        end
    end
    input = opts["input"]
    output = opts["output"]
    isempty(input) && error("--input is required")
    abspath(input) == abspath(output) &&
        error("--output must differ from --input; both resolve to $(abspath(input))")
    dsd = opts["dsd"]
    dsd in ("single", "double") || error("--dsd must be single or double")
    (; input,
        output,
        dsd,
        freq_hz = parse(Float64, opts["freq"]),
        diameters = parse_range(opts["diameters"]),
        cache = opts["cache"],
        solver = IITM(parse(Int, opts["nmax"]), parse(Int, opts["nr"]),
            parse(Int, opts["ntheta"])))
end

function parse_range(spec)
    parts = parse.(Float64, split(spec, ":"))
    bins = if length(parts) == 3
        parts[2] != 0 || error("range step must be nonzero")
        collect(parts[1]:parts[2]:parts[3])
    elseif length(parts) == 2
        collect(parts[1]:parts[2])
    else
        error("range must be start:step:stop or start:stop")
    end
    isempty(bins) && error("range must produce at least one bin")
    bins
end

function water_refractive_index(freq_hz, T_celsius)
    λ_m = C0_M_PER_S / freq_hz
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
            σ * λ_m / (2π * C0_M_PER_S * ε0)

    sqrt(complex(εreal, εimag))
end

function rain_axis_ratio(D_mm)
    D_mm < 0.7 && return 1.0
    D_mm < 1.5 && return 1.173 - 0.5265D_mm + 0.4698D_mm^2 -
                        0.1317D_mm^3 - 8.5e-3D_mm^4
    1.065 - 6.25e-2D_mm - 3.99e-3D_mm^2 +
        7.66e-4D_mm^3 - 4.095e-5D_mm^4
end

function scattering_table(Ds, freq_hz; temperature_celsius = 0.0,
        solver = IITM(6, 10, 16))
    λ = C0_MM_PER_S / freq_hz
    m = water_refractive_index(freq_hz, temperature_celsius)
    map(Ds) do D
        a = D / 2
        c = a * rain_axis_ratio(D)
        shape = Spheroid{Float64, ComplexF64}(a, c, ComplexF64(m))
        T = transition_matrix(shape, λ, solver)
        Sback = amplitude_matrix(T, π / 2, 0.0, π / 2, π; λ)
        Sfwd = amplitude_matrix(T, π / 2, 0.0, π / 2, 0.0; λ)
        (; D_mm = D,
            Shh_back = Sback[2, 2],
            Svv_back = Sback[1, 1],
            Shh_forward = Sfwd[2, 2],
            Svv_forward = Sfwd[1, 1])
    end
end

scattering_solver_params(solver) = (; solver = repr(solver))

function scattering_solver_params(solver::IITM)
    (; nmax = solver.nₘₐₓ,
        nr = solver.Nr,
        ntheta = solver.Nϑ,
        nphi = solver.Nφ,
        rmin = solver.rₘᵢₙ)
end

function scattering_cache_matches(payload, Ds, freq_hz, solver_params)
    names = cache_propertynames(payload)
    (:Ds in names && :freq_hz in names && :solver_params in names) || return false
    payload.Ds == Ds && payload.freq_hz == freq_hz &&
        payload.solver_params == solver_params
end

function cache_propertynames(payload)
    try
        propertynames(payload)
    catch
        ()
    end
end

function cached_scattering_table(payload, Ds)
    :table in cache_propertynames(payload) || return nothing
    table = payload.table
    isnothing(table) && return nothing
    try
        length(table) == length(Ds) ? table : nothing
    catch
        nothing
    end
end

function load_or_compute_scattering_table(cache, Ds, freq_hz; solver)
    solver_params = scattering_solver_params(solver)
    if !isempty(cache) && isfile(cache)
        payload = try
            deserialize(cache)
        catch err
            @warn "ignoring unreadable scattering table cache" cache error = sprint(showerror, err)
            nothing
        end
        if !isnothing(payload)
            if scattering_cache_matches(payload, Ds, freq_hz, solver_params)
                table = cached_scattering_table(payload, Ds)
                if !isnothing(table)
                    println("loaded scattering table cache: ", cache)
                    return table
                end
                @warn "ignoring scattering table cache without a valid table" cache
            else
                println("ignoring stale scattering table cache: ", cache)
            end
        end
    end

    table = scattering_table(Ds, freq_hz; solver)
    if !isempty(cache)
        serialize(cache, (; Ds, freq_hz, solver_params, table))
        println("wrote scattering table cache: ", cache)
    end
    table
end

struct WRFFields
    xlat::Array{Float64, 3}
    xlong::Array{Float64, 3}
    znu::Vector{Float64}
    xtime::Vector{Float64}
    qrain::Array{Float64, 4}
    qnrain::Union{Nothing, Array{Float64, 4}}
    qvapor::Array{Float64, 4}
    pb::Array{Float64, 4}
    p::Array{Float64, 4}
    theta::Array{Float64, 4}
end

read_array(ds, name) = Float64.(ds[name][ntuple(_ -> Colon(), ndims(ds[name]))...])

function read_wrf_fields(path; dsd)
    ds = NCDataset(path, "r")
    try
        qnrain = haskey(ds, "QNRAIN") ? read_array(ds, "QNRAIN") : nothing
        dsd == "double" && isnothing(qnrain) &&
            error("double DSD mode requires QNRAIN in the input file")
        WRFFields(
            read_array(ds, "XLAT"),
            read_array(ds, "XLONG"),
            vec(read_array(ds, "ZNU")),
            vec(read_array(ds, "XTIME")),
            read_array(ds, "QRAIN"),
            qnrain,
            read_array(ds, "QVAPOR"),
            read_array(ds, "PB"),
            read_array(ds, "P"),
            read_array(ds, "T"))
    finally
        close(ds)
    end
end

function dsd_parameters(qrv, qn_density, mode)
    if mode == "single"
        Λ = (MASS_COEFF * RAIN_N0 * gamma(4) / qrv)^(1 / 4)
        return (; N0 = RAIN_N0, Λ)
    end

    Λ = (MASS_COEFF * gamma(4) * qn_density / qrv)^(1 / 3)
    N0 = qn_density * Λ
    (; N0, Λ)
end

function radar_moments(table, freq_hz, N0, Λ)
    λ = C0_MM_PER_S / freq_hz
    ΔD = length(table) > 1 ? table[2].D_mm - table[1].D_mm : 1.0
    m = water_refractive_index(freq_hz, 0.0)
    Kw2 = abs2((m^2 - 1) / (m^2 + 2))
    weights = [N0 * exp(-Λ * row.D_mm) * ΔD for row in table]

    hh = sum(abs2(row.Shh_back) * weights[i] for (i, row) in pairs(table))
    vv = sum(abs2(row.Svv_back) * weights[i] for (i, row) in pairs(table))
    if !(hh > 0 && vv > 0)
        return (; ZH = NaN, ZV = NaN, ZDR = NaN, RHOHV = NaN, KDP = NaN)
    end

    Zhh = hh * λ^4 / (π^5 * Kw2)
    Zvv = vv * λ^4 / (π^5 * Kw2)
    cov = sum(conj(row.Shh_back) * row.Svv_back * weights[i]
        for (i, row) in pairs(table))
    ρhv = abs(cov) / sqrt(hh * vv)
    kdp = 180 / π * 1e-3 * λ *
          sum(real(row.Shh_forward - row.Svv_forward) * weights[i]
              for (i, row) in pairs(table))

    (; ZH = 10log10(Zhh),
        ZV = 10log10(Zvv),
        ZDR = 10log10(Zhh / Zvv),
        RHOHV = ρhv,
        KDP = kdp)
end

function valid_radar_moments(moments)
    all(isfinite, (moments.ZH, moments.ZV, moments.ZDR, moments.RHOHV, moments.KDP))
end

function process_wrf_grid(fields::WRFFields, table, freq_hz; dsd)
    shape = size(fields.qrain)
    ZH = fill(FILL, shape)
    ZV = fill(FILL, shape)
    ZDR = fill(FILL, shape)
    RHOHV = fill(FILL, shape)
    KDP = fill(FILL, shape)
    QNRAIN_OUT = isnothing(fields.qnrain) ? fill(FILL, shape) : Float32.(fields.qnrain)

    Threads.@threads for idx in CartesianIndices(fields.qrain)
        qr = fields.qrain[idx]
        qr > 1e-10 || continue

        rho = air_density(fields.qvapor[idx], fields.pb[idx], fields.p[idx],
            fields.theta[idx])
        qrv = qr * rho
        qn_density = dsd == "double" ? max(fields.qnrain[idx] * rho, 1e-9) : NaN
        params = dsd_parameters(qrv, qn_density, dsd)
        dsd == "single" && (QNRAIN_OUT[idx] = Float32(params.N0 / params.Λ / rho))

        moments = radar_moments(table, freq_hz, params.N0, params.Λ)
        if valid_radar_moments(moments)
            ZH[idx] = Float32(moments.ZH)
            ZV[idx] = Float32(moments.ZV)
            ZDR[idx] = Float32(moments.ZDR)
            RHOHV[idx] = Float32(moments.RHOHV)
            KDP[idx] = Float32(moments.KDP)
        else
            ZH[idx] = FILL
            ZV[idx] = FILL
            ZDR[idx] = FILL
            RHOHV[idx] = FILL
            KDP[idx] = FILL
        end
    end

    (; ZH, ZV, ZDR, RHOHV, KDP, QNRAIN = QNRAIN_OUT)
end

function write_output(path, fields::WRFFields, output)
    rm(path; force = true)
    ds = NCDataset(path, "c")
    try
        nx, ny, nz, nt = size(fields.qrain)
        defDim(ds, "west_east", nx)
        defDim(ds, "south_north", ny)
        defDim(ds, "bottom_top", nz)
        defDim(ds, "Time", nt)
        ds.attrib["title"] = "Synthetic dual-pol radar output"
        ds.attrib["source"] = "TransitionMatrices.jl radiative_transfer example"

        write_var(ds, "XLAT", Float32.(fields.xlat), ("west_east", "south_north", "Time"))
        write_var(ds, "XLONG", Float32.(fields.xlong), ("west_east", "south_north", "Time"))
        write_var(ds, "ZNU", Float32.(fields.znu), ("bottom_top",))
        write_var(ds, "XTIME", Float32.(fields.xtime), ("Time",))

        dims4 = ("west_east", "south_north", "bottom_top", "Time")
        write_var(ds, "ZH", output.ZH, dims4; units = "dBZ")
        write_var(ds, "ZV", output.ZV, dims4; units = "dBZ")
        write_var(ds, "ZDR", output.ZDR, dims4; units = "dB")
        write_var(ds, "RHOHV", output.RHOHV, dims4)
        write_var(ds, "KDP", output.KDP, dims4; units = "degree km-1")
        write_var(ds, "QRAIN", Float32.(fields.qrain), dims4; units = "kg kg-1")
        write_var(ds, "QNRAIN", output.QNRAIN, dims4; units = "kg-1")
    finally
        close(ds)
    end
end

function main(args = ARGS)
    opts = parse_args(args)
    fields = read_wrf_fields(opts.input; dsd = opts.dsd)
    table = load_or_compute_scattering_table(opts.cache, opts.diameters, opts.freq_hz;
        solver = opts.solver)
    output = process_wrf_grid(fields, table, opts.freq_hz; dsd = opts.dsd)
    write_output(opts.output, fields, output)
    valid = count(!=(FILL), output.ZH)
    println("wrote dual-pol radar output: ", opts.output)
    println("valid rain grid cells: ", valid, " / ", length(output.ZH))
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
