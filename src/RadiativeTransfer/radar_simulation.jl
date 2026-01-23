# This simulator follows the methodology of:
# - H. Li, Y. Xiong and Y. Chen, "Simulation of Complex Meteorological
#   Target Echoes for Airborne Dual-Polarization Weather Radar Based
#   on Invariant Imbedding T -Matrix," in IEEE Transactions on Geoscience
#   and Remote Sensing, vol. 62, pp. 1–17, 2024, Art. no. 5105817.
# "Validation of Simulated Hurricane Drop Size Distributions
# using Polarimetric Radar", Geophys. Res. Lett., 42,
# doi:10.1002/2015GL067278.
#Download WRF sample data from:  
#🔗 **Link:** https://pan.baidu.com/s/1xFomb3d1TYG66LD_NZMghw  
#**Extraction code:** `hbub`

using SpecialFunctions
using ArgParse
using NetCDF
using DataStructures
using TransitionMatrices
using LinearAlgebra
using Printf
using StaticArrays
using MAT  # For saving/loading .mat files

# Global constants
xam_r = pi * 997.0 / 6.0
xbm_r = 3.0
xmu_r = 0.0
xobmr = 1.0 / xbm_r
xcre = 1.0 + 2 * xbm_r + xmu_r

# Function to parse command line arguments
function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--dsdtype", "-d"
        help = "single, double, full"
        arg_type = String
        required = true
        "--output", "-o"
        help = "path to output file"
        arg_type = String
        default = "dualpol_radar.nc"
        "input"
        help = "path to input file"
        arg_type = String
        required = true
    end
    return parse_args(s)
end

# Parameters for the calculation
D_vals = range(0.5, 8.5, length=50)  # Diameter range in μm
nmax = 30                            # Maximum order for Mie calculation
k = 2 * π /32                            # Wavenumber, assuming λ = 1
λ = 32
m = 7.253701208520992 + 2.855178298595981im  # Complex permittivity for calculation

# Function to calculate radar variables once and save them
function calculate_and_save_scattering(D_vals, nmax, m, λ, filename)
    # Lists to store scattering parameters
    S11_list_b = ComplexF64[]  # HH for backward scattering
    S12_list_b = ComplexF64[]  # HV for backward scattering
    S21_list_b = ComplexF64[]  # VH for backward scattering
    S22_list_b = ComplexF64[]  # VV for backward scattering

    S11_list_f = ComplexF64[]  # HH for forward scattering
    S12_list_f = ComplexF64[]  # HV for forward scattering
    S21_list_f = ComplexF64[]  # VH for forward scattering
    S22_list_f = ComplexF64[]  # VV for forward scattering
    function axis_ratio(D::Float64)
        if D < 0.7
            return 1.0
        elseif D >= 0.7 && D < 1.5
            return 1.173 - 0.5265 * D + 0.4698 * D^2 - 0.1317 * D^3 - 8.5e-3 * D^4
        else
            return 1.065 - 6.25e-2 * D - 3.99e-3 * D^2 + 7.66e-4 * D^3 - 4.095e-5 * D^4
        end
    end
    # Main loop over particle diameter
    for D in D_vals
        r = D / 2
        # ar = axis_ratio(D)
        ar = 1
        s = Spheroid{Float64, ComplexF64}(r,ar*r,m)
        Tm =calc_T(s,λ)

        # Backward scattering (φ = π)
        S_b = calc_S(Tm, 1.57, 0, 1.57, π; λ = λ)
        push!(S11_list_b, S_b[1, 1])
        push!(S12_list_b, S_b[1, 2])
        push!(S21_list_b, S_b[2, 1])
        push!(S22_list_b, S_b[2, 2])

        # Forward scattering (φ = 0)
        S_f = calc_S(Tm, 1.57, 0, 1.57, 0; λ = λ)
        push!(S11_list_f, S_f[1, 1])
        push!(S12_list_f, S_f[1, 2])
        push!(S21_list_f, S_f[2, 1])
        push!(S22_list_f, S_f[2, 2])
    end

    # Save scattering matrices to .mat file for later use
    scattering_data = Dict(
        "S11_b" => S11_list_b,
        "S12_b" => S12_list_b,
        "S21_b" => S21_list_b,
        "S22_b" => S22_list_b,
        "S11_f" => S11_list_f,
        "S12_f" => S12_list_f,
        "S21_f" => S21_list_f,
        "S22_f" => S22_list_f
    )

    matwrite(filename, scattering_data)
    println("Scattering matrices saved to ", filename)
end

# Function to load scattering matrices from a .mat file
function load_scattering_matrices(filename)
    scattering_data = matread(filename)
    return scattering_data["S11_b"], scattering_data["S12_b"], scattering_data["S21_b"], scattering_data["S22_b"],
           scattering_data["S11_f"], scattering_data["S12_f"], scattering_data["S21_f"], scattering_data["S22_f"]
end

# Function to read variables from WRF NetCDF file
function read_nc_var(filename, varnames::Array)
    data = OrderedDict()
    for varname in varnames
        data[string(varname)] = ncread(filename, string(varname))
    end
    return collect(values(data))
end

# Function to save data to MAT file
function write_matfile(filename, lat, lon, lev, times, ZH, ZDR, ZV, ph,Kdp, QRAIN, QNRAIN)
    println("Write to mat file ...")
    data_dict = Dict(
        "XLAT" => lat,
        "XLON" => lon,
        "ZNU" => lev,
        "XTIME" => times,
        "ZH" => ZH,
        "ZDR" => ZDR,
        "ZV" => ZV,
        "ph" => ph,
        "Kdp" => Kdp,
        "QRAIN" => QRAIN,
        "QNRAIN" => QNRAIN
    )
    matwrite(filename, data_dict)
    println("Saved to ", filename)
end

# Function to calculate radar variables using the prescribed gamma distribution
function calc_radar_gamma_dsd(N0, lambda)
    m = 7.253701208520992 + 2.855178298595981im

    λ = 32  # Wavelength in meters
    Kw = (m - 1) / (m+ 2)  # Dielectric factor
    Kw2 = abs2(Kw)
    ΔD = step(D_vals)

    Λ = lambda        # Slope parameter (mm⁻¹)
    N_D = D -> N0 * exp(-Λ * D)

    σ_hh_list = abs2.(S22_list_b)
    σ_vv_list = abs2.(S11_list_b)

    Zhh = sum(σ_hh_list[i] * N_D(D_vals[i]) * ΔD for i in eachindex(D_vals))
    Zvv = sum(σ_vv_list[i] * N_D(D_vals[i]) * ΔD for i in eachindex(D_vals))

    # Apply radar reflectivity constant
    C = λ^4 / (π^5 * Kw2)
    Zhh *= C
    Zvv *= C
    Zdr = 10 * log10(Zhh / Zvv)

    # Convert Zhh and Zvv to dBZ
    dBZhh = 10 * log10(Zhh)
    dBZvv = 10 * log10(Zvv)


    # ----------------- Computation of ρ_h (correlation coefficient) -----------------

    numerator = 0.0 + 0.0im
    denom1 = 0.0
    denom2 = 0.0

    for i in eachindex(D_vals)
        nD = N_D(D_vals[i])
        Shh = S11_list_b[i]
        Svv = S22_list_b[i]

        numerator += conj(Shh) * Svv * nD * ΔD
        denom1 += abs2(Shh) * nD * ΔD
        denom2 += abs2(Svv) * nD * ΔD
    end

    ρh = abs(numerator) / sqrt(denom1 * denom2)


    # ----------------- Computation of Kdp (differential phase shift rate) -----------------

    Kdp_integral = 0.0
    for i in eachindex(D_vals)
        Shh = S22_list_f[i]
        Svv = S11_list_f[i]
        diff_real = real(Shh - Svv)
        Kdp_integral += diff_real * N_D(D_vals[i]) * ΔD
    end

    Kdp = (180 / π) *1e-3* λ * Kdp_integral
    # println(Kdp)
    return dBZhh,dBZvv,dBZhh-dBZvv,ρh,Kdp

end

# Main program
args = parse_commandline()
filein = args["input"]
fileout = args["output"]
dsdtype = args["dsdtype"]

# Check if the scattering matrices have already been computed
scattering_file = "scattering_matrices.mat"  # File to store the scattering matrices

if isfile(scattering_file)
    println("Loading precomputed scattering matrices from ", scattering_file)
    S11_list_b, S12_list_b, S21_list_b, S22_list_b,
    S11_list_f, S12_list_f, S21_list_f, S22_list_f = load_scattering_matrices(scattering_file)
else
    println("Calculating and saving scattering matrices...")
    calculate_and_save_scattering(D_vals, nmax, m, λ, scattering_file)
    # Load the newly computed scattering matrices
    S11_list_b, S12_list_b, S21_list_b, S22_list_b,
    S11_list_f, S12_list_f, S21_list_f, S22_list_f = load_scattering_matrices(scattering_file)
end

# Read the WRF data (assuming variables like lat, lon, lev, times, etc.)
println("Reading WRF data from ", filein)
fileinfo = ncinfo(filein)

dims = ["XLAT", "XLONG", "ZNU", "XTIME"]
lat, lon, lev, times = read_nc_var(filein, dims)

if dsdtype == "single"
    vars = ["QRAIN", "QVAPOR", "PB", "P", "T", "REFL_10CM"]
    qrain, qv, pb, pp, theta, refl = read_nc_var(filein, vars)
elseif dsdtype == "double"
    vars = ["QRAIN", "QNRAIN", "QVAPOR", "PB", "P", "T", "REFL_10CM"]
    qrain, qnrain, qv, pb, pp, theta, refl = read_nc_var(filein, vars)
elseif dsdtype == "full"
    vars = ["QRAIN", "QNRAIN", "QVAPOR", "PB", "P", "T", "REFL_10CM", "ff1i01", "ff1i02", "ff1i03", "ff1i04", "ff1i05"]
    ff1 = Vector{Array{Float32, 4}}(undef, 33)
    for d in 1:33
        ff1[d] = ncread(filein, "ff1i$(d)")
    end
    qrain, qnrain, qv, pb, pp, theta, refl = read_nc_var(filein, vars)
else
    println(dsdtype, " not recognized")
    exit
end

# Initialize radar variables
ZDR = similar(qrain); fill!(ZDR, -999.0)
ZV = similar(qrain); fill!(ZV, -999.0)
ZH = similar(qrain); fill!(ZH, -999.0)
ρh = similar(qrain); fill!(ρh, -999.0)
Kdp = similar(qrain); fill!(Kdp, -999.0)

# Calculate radar variables for each WRF grid cell
println("Calculating radar variables for each WRF grid cell...")
gamma1 = gamma(1. + xbm_r + xmu_r)
gamma2 = gamma(1. + xmu_r)
xorg2 = 1.0 / gamma2

Threads.@threads for i in 1:size(qrain, 1)  # Parallelizing over the i dimension
    for j in 1:size(qrain, 2)
        for k in 1:size(qrain, 3)
            for t in 1:size(qrain, 4)  # No need to parallelize this loop as it is only 1 in your case
                qvapor = max(1.0e-10, qv[i, j, k, t])
                pressure = pb[i, j, k, t] + pp[i, j, k, t]
                tempk = (theta[i, j, k, t] + 300.0) * (pressure / 100000.0)^(2.0 / 7.0)
                rho = 0.622 * pressure / (287.15 * tempk * (qvapor + 0.622))

                if (qrain[i, j, k, t] > 1.0e-9)
                    qrv = qrain[i, j, k, t] * rho
                    N0 = 8000
                    lambda = ((xam_r * gamma1 * N0 * 1000.0 / qrv)^(1.0 / (1.0 + xbm_r))) / 1000.0
                    qnrain[i, j, k, t] = N0 * gamma2 / (lambda^(1.0 + xmu_r))
                    # Call the radar calculation function for ZH, ZV, ZDR, ρh, Kdp
                    ZH[i, j, k, t], ZV[i, j, k, t], ZDR[i, j, k, t], ρh[i, j, k, t], Kdp[i, j, k, t] = calc_radar_gamma_dsd(N0, lambda)
                else

                end
            end
        end
    end
end

# Write the results to a .mat file
write_matfile(fileout, lat, lon, lev, times, ZH, ZDR, ZV, ρh, Kdp, qrain, qnrain)
