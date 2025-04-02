# refractive_index.jl
#
# This module provides functions to compute the complex refractive index
# of water and ice as a function of frequency and temperature.
#
# The refractive index model is primarily based on the Debye-type formulation
# for liquid water and standard models for ice in microwave frequencies.
#
# The implementation follows the formulation described in:
# G. Zhang, *Weather Radar Polarimetry*, Beijing, China: China Meteorological Press, 2018, pp. 39–40.
#
# These models are commonly used in microwave remote sensing and radar polarimetry
# simulations for atmospheric hydrometeors.
#


# Define constants
T = 0.0
seita_0 = 8.854e-12
c = 3e8
freq = 9.375e9
lamda = c / freq

# Water parameters (Debye formula)
kesei_s = 78.54 * (1.0 - 4.579e-3 * (T - 25) + 1.19e-5 * ((T - 25)^2) - 2.8e-8 * ((T - 25))^3)
kesei_oo = 5.27137 + 2.16474e-2 * T - 1.31198e-3 * (T^2)
alpha = (-16.8129 / (T + 273)) + 6.09265e-2
lamda_s = 3.3836e-6 * exp(2513.98 / (T + 273))
seita = 1.1117e-4

kesei_real = kesei_oo + ((kesei_s - kesei_oo) * (1 + (lamda_s / lamda)^(1 - alpha) * sin(alpha * π / 2))) / 
                (1 + 2 * ((lamda_s / lamda)^(1 - alpha) * sin(alpha * π / 2)) + (lamda_s / lamda)^(2 - 2 * alpha))
kesei_imag = (((kesei_s - kesei_oo) * ((lamda_s / lamda)^(1 - alpha) * cos(alpha * π / 2))) /
                (1 + 2 * ((lamda_s / lamda)^(1 - alpha) * sin(alpha * π / 2)) + (lamda_s / lamda)^(2 - 2 * alpha))) +
                (seita * lamda) / (2 * π * c * seita_0)

W = sqrt(kesei_real + kesei_imag * im)
println(W)
k_water = W * W

# Ice parameters
T = 0.0
kesei_s_ice = 203.168 + 2.5 * T + 0.15 * T^2
kesei_oo_ice = 3.168
alpha_ice = 0.288 + 5.2e-3 * T + 2.3e-4 * T^2
lamda_s_ice = 9.990288e-6 * exp(6643.5 / (T + 273))
seita_ice = 1.1146e-13 * exp(-6291.2 / (T + 273))

kesei_real_ice = kesei_oo_ice + ((kesei_s_ice - kesei_oo_ice) * (1 + (lamda_s_ice / lamda)^(1 - alpha_ice) * sin(alpha_ice * π / 2))) / 
                        (1 + 2 * ((lamda_s_ice / lamda)^(1 - alpha_ice) * sin(alpha_ice * π / 2)) + (lamda_s_ice / lamda)^(2 - 2 * alpha_ice))

kesei_imag_ice = (((kesei_s_ice - kesei_oo_ice) * ((lamda_s_ice / lamda)^(1 - alpha_ice) * cos(alpha_ice * π / 2))) /
                    (1 + 2 * ((lamda_s_ice / lamda)^(1 - alpha_ice) * sin(alpha_ice * π / 2)) + (lamda_s_ice / lamda)^(2 - 2 * alpha_ice))) +
                    (seita_ice * lamda) / (2 * π * c * seita_0)

W_ice = sqrt(kesei_real_ice + kesei_imag_ice * im)
k_ice = W_ice * W_ice
