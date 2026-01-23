using TransitionMatrices
using MAT
using Printf
using MAT
using Printf
# Constants
c = 3e8
T = 0.0
D_list = collect(0.1:0.2:10.0)       # Diameter list in mm
frequencies = [10e9, 18e9]          # Frequencies in Hz

# Complex refractive index of water (Debye model)
function water_refractive_index(freq::Float64, T::Float64)
    seita_0 = 8.854e-12
    c = 3e8
    lamda = c / freq
    kesei_s = 78.54 * (1.0 - 4.579e-3 * (T - 25) + 1.19e-5 * ((T - 25)^2) - 2.8e-8 * ((T - 25)^3))
    kesei_oo = 5.27137 + 2.16474e-2 * T - 1.31198e-3 * (T^2)
    alpha = (-16.8129 / (T + 273)) + 6.09265e-2
    lamda_s = 3.3836e-6 * exp(2513.98 / (T + 273))
    seita = 1.1117e-4
    num = (lamda_s / lamda)^(1 - alpha)
    denom = 1 + 2*num*sin(alpha*pi/2) + (lamda_s / lamda)^(2 - 2*alpha)
    kesei_real = kesei_oo + ((kesei_s - kesei_oo) * (1 + num * sin(alpha*pi/2))) / denom
    kesei_imag = ((kesei_s - kesei_oo) * (num * cos(alpha*pi/2))) / denom + (seita * lamda) / (2 * π * c * seita_0)
    return sqrt(complex(kesei_real, kesei_imag))
end

# Initialize extinction coefficient lookup table
Cext_table = zeros(length(frequencies), length(D_list))

# Main loop for computing extinction cross-sections
for (j, f) in enumerate(frequencies)
    λ = c / f * 1000
    m = water_refractive_index(f, T)
    @printf("Freq %.1f GHz, m = %.4f + %.4fi\n", f/1e9, real(m), imag(m))

    for (i, D) in enumerate(D_list)
        a = D / 2                     # Semi-major axis
        b = a * (1 - 0.062 * D)       # Semi-minor axis
        s = Spheroid(a, b, m)
        try
            Tmat = calc_T(s, λ)
            Cext_table[j, i] = calc_Cext(Tmat)
        catch
            Cext_table[j, i] = 0.0
        end
    end
end

# Save extinction cross-section table
matwrite("Cext_table_rain.mat", Dict(
    "frequencies" => frequencies,
    "diameters" => D_list,
    "Cext_table" => Cext_table
))
println("✅ Saved Cext_table_rain.mat")

# ---------------------------------------------------------------------

# Parameters
T_atm = 250.0           # Atmospheric temperature in K
H = 10.0                # Path length in km
N0 = 8e3                # Intercept parameter (m⁻³ mm⁻¹)
rain_rates = collect(0:10:100)  # Rainfall rates in mm/h
ϵ = 1e-6

# Load extinction data
data = matread("Cext_table_rain.mat")
frequencies = data["frequencies"]
D_list = data["diameters"]
Cext_table = data["Cext_table"]
dr = D_list[2] - D_list[1]

# Initialize brightness temperature (TB) lookup table
TB_LUT = zeros(length(rain_rates), length(frequencies))

# Compute TB values for each rain rate and frequency
for (j, f) in enumerate(frequencies)
    for (i, R) in enumerate(rain_rates)
        if R == 0
            TB_LUT[i, j] = 0.0
            continue
        end
        Λ = 4.1 * R^(-0.21)
        Nr = N0 .* exp.(-Λ .* D_list)                 # Drop size distribution (m⁻³ mm⁻¹)
        α_ext = sum(Nr .* Cext_table[j, :]) * dr / 1e6  # Convert to km⁻¹
        println(α_ext)
        τ = α_ext * H                                  # Optical depth
        TB_LUT[i, j] = T_atm * (1 - exp(-τ))           # Radiative transfer equation
        # @printf "R = %.1f mm/h, freq = %.1f GHz → α = %.4f, TB = %.2f K\n" R, f/1e9, α_ext, TB_LUT[i, j]
    end
end

# Save brightness temperature lookup table
matwrite("TB_LUT_Rain_From_Cext.mat", Dict(
    "frequencies" => frequencies,
    "rain_rates" => rain_rates,
    "TB_LUT" => TB_LUT
))
println("✅ Brightness temperature LUT saved to TB_LUT_Rain_From_Cext.mat")
