#!/usr/bin/env julia

using NCDatasets

include(joinpath(@__DIR__, "shared.jl"))
using .RadiativeTransferShared: MASS_COEFF, RAIN_N0, R_DRY_AIR, WATER_DENSITY,
                                air_density, write_var

function usage()
    println("""
    Generate a tiny deterministic WRF-like NetCDF file for smoke testing.

    Usage:
      julia --project=examples/radiative_transfer \\
        examples/radiative_transfer/generate_fake_wrf.jl [options]

    Options:
      --output PATH   Output NetCDF path [fake_wrfout.nc]
      --nx N          west_east grid size [8]
      --ny N          south_north grid size [6]
      --nz N          bottom_top grid size [4]
      --nt N          time size [1]
      --help          Show this message
    """)
end

function parse_args(args)
    opts = Dict(
        "output" => "fake_wrfout.nc",
        "nx" => "8",
        "ny" => "6",
        "nz" => "4",
        "nt" => "1")
    i = 1
    while i <= length(args)
        arg = args[i]
        if arg == "--help" || arg == "-h"
            usage()
            exit(0)
        elseif haskey(opts, replace(arg, "--" => ""))
            i == length(args) && error("missing value for $arg")
            opts[replace(arg, "--" => "")] = args[i + 1]
            i += 2
        else
            error("unknown argument: $arg")
        end
    end
    nx = parse(Int, opts["nx"])
    ny = parse(Int, opts["ny"])
    nz = parse(Int, opts["nz"])
    nt = parse(Int, opts["nt"])
    for (name, value) in (("nx", nx), ("ny", ny), ("nz", nz), ("nt", nt))
        value > 0 || error("invalid --$name: expected a positive integer, got $value")
    end
    (; output = opts["output"], nx, ny, nz, nt)
end

function synthetic_fields(nx, ny, nz, nt)
    xlat = Array{Float32}(undef, nx, ny, nt)
    xlon = similar(xlat)
    znu = collect(Float32, range(1.0, 0.15; length = nz))
    xtime = collect(Float32, 0:10:(10 * (nt - 1)))

    shape = (nx, ny, nz, nt)
    qrain = zeros(Float32, shape)
    qnrain = zeros(Float32, shape)
    qvapor = Array{Float32}(undef, shape)
    pb = Array{Float32}(undef, shape)
    p = Array{Float32}(undef, shape)
    theta = Array{Float32}(undef, shape)
    refl = fill(Float32(-999), shape)

    cx = (nx + 1) / 2
    cy = (ny + 1) / 2
    sx = max(nx / 4, 1)
    sy = max(ny / 4, 1)

    for t in 1:nt
        for j in 1:ny, i in 1:nx

            xlat[i, j, t] = 34.0f0 + 0.04f0 * Float32(j - 1)
            xlon[i, j, t] = -98.0f0 + 0.04f0 * Float32(i - 1)
        end

        for k in 1:nz, j in 1:ny, i in 1:nx
            height_factor = exp(-0.45 * (k - 1))
            blob = exp(-((i - cx)^2 / (2sx^2) + (j - cy)^2 / (2sy^2) +
                         (k - 1)^2 / 2.5))
            rain = blob > 0.035 ? 8.0e-4 * blob * height_factor : 0.0

            qvapor[i, j, k, t] = Float32(0.011 * exp(-0.12 * (k - 1)))
            pb[i, j, k, t] = Float32(93000 * exp(-0.10 * (k - 1)))
            p[i, j, k, t] = Float32(250 * sinpi((i - 1) / max(nx - 1, 1)) *
                                    cospi((j - 1) / max(ny - 1, 1)))
            theta[i, j, k, t] = Float32(2.0 - 1.5 * (k - 1))
            qrain[i, j, k, t] = Float32(rain)
            refl[i, j, k, t] = Float32(rain > 0 ? 10 + 45 * sqrt(blob) : -999)

            if rain > 0
                rho = air_density(qvapor[i, j, k, t], pb[i, j, k, t],
                    p[i, j, k, t], theta[i, j, k, t])
                qrv = rain * rho
                lambda = (MASS_COEFF * RAIN_N0 * 6 / qrv)^(1 / 4)
                number_density = RAIN_N0 / lambda
                qnrain[i, j, k, t] = Float32(number_density / rho)
            end
        end
    end

    (; xlat, xlon, znu, xtime, qrain, qnrain, qvapor, pb, p, theta, refl)
end

function write_fake_wrf(path, fields)
    rm(path; force = true)
    ds = NCDataset(path, "c")
    try
        nx, ny, nz, nt = size(fields.qrain)
        defDim(ds, "west_east", nx)
        defDim(ds, "south_north", ny)
        defDim(ds, "bottom_top", nz)
        defDim(ds, "Time", nt)
        ds.attrib["title"] = "Synthetic WRF-like smoke-test data"
        ds.attrib["source"] = "TransitionMatrices.jl examples/radiative_transfer"
        ds.attrib["note"] = "Not validation data; generated only to exercise the pipeline."

        write_var(ds, "XLAT", fields.xlat, ("west_east", "south_north", "Time");
            units = "degree_north")
        write_var(ds, "XLONG", fields.xlon, ("west_east", "south_north", "Time");
            units = "degree_east")
        write_var(ds, "ZNU", fields.znu, ("bottom_top",);
            description = "synthetic eta coordinate")
        write_var(ds, "XTIME", fields.xtime, ("Time",); units = "minutes")

        dims4 = ("west_east", "south_north", "bottom_top", "Time")
        write_var(ds, "QRAIN", fields.qrain, dims4; units = "kg kg-1")
        write_var(ds, "QNRAIN", fields.qnrain, dims4; units = "kg-1")
        write_var(ds, "QVAPOR", fields.qvapor, dims4; units = "kg kg-1")
        write_var(ds, "PB", fields.pb, dims4; units = "Pa")
        write_var(ds, "P", fields.p, dims4; units = "Pa")
        write_var(ds, "T", fields.theta, dims4; units = "K",
            description = "WRF-style perturbation potential temperature")
        write_var(ds, "REFL_10CM", fields.refl, dims4; units = "dBZ",
            description = "synthetic placeholder reflectivity")
    finally
        close(ds)
    end
end

function main(args = ARGS)
    opts = parse_args(args)
    fields = synthetic_fields(opts.nx, opts.ny, opts.nz, opts.nt)
    write_fake_wrf(opts.output, fields)
    println("wrote synthetic WRF-like file: ", opts.output)
    println("grid: $(opts.nx) x $(opts.ny) x $(opts.nz) x $(opts.nt)")
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
