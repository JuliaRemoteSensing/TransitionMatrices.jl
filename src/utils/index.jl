"""
`phase_matrix(𝐒::AbstractMatrix)`

Calculate the phase matrix `𝐙` from the amplitude matrix `𝐒`, according to Eq. (2.106) -- Eq. (2.121) in Mishchenko et al. (2002).
"""
function phase_matrix(𝐒::AbstractMatrix)
    𝐙₁₁ = 0.5 * (𝐒[1, 1] * 𝐒[1, 1]' + 𝐒[1, 2] * 𝐒[1, 2]' + 𝐒[2, 1] * 𝐒[2, 1]' +
           𝐒[2, 2] * 𝐒[2, 2]')
    𝐙₁₂ = 0.5 * (𝐒[1, 1] * 𝐒[1, 1]' - 𝐒[1, 2] * 𝐒[1, 2]' + 𝐒[2, 1] * 𝐒[2, 1]' -
           𝐒[2, 2] * 𝐒[2, 2]')
    𝐙₁₃ = -𝐒[1, 1] * 𝐒[1, 2]' - 𝐒[2, 2] * 𝐒[2, 1]'
    𝐙₁₄ = 1.0im * (𝐒[1, 1] * 𝐒[1, 2]' - 𝐒[2, 2] * 𝐒[2, 1]')
    𝐙₂₁ = 0.5 * (𝐒[1, 1] * 𝐒[1, 1]' + 𝐒[1, 2] * 𝐒[1, 2]' - 𝐒[2, 1] * 𝐒[2, 1]' -
           𝐒[2, 2] * 𝐒[2, 2]')
    𝐙₂₂ = 0.5 * (𝐒[1, 1] * 𝐒[1, 1]' - 𝐒[1, 2] * 𝐒[1, 2]' - 𝐒[2, 1] * 𝐒[2, 1]' +
           𝐒[2, 2] * 𝐒[2, 2]')
    𝐙₂₃ = -𝐒[1, 1] * 𝐒[1, 2]' + 𝐒[2, 2] * 𝐒[2, 1]'
    𝐙₂₄ = 1.0im * (𝐒[1, 1] * 𝐒[1, 2]' + 𝐒[2, 2] * 𝐒[2, 1]')
    𝐙₃₁ = -𝐒[1, 1] * 𝐒[2, 1]' - 𝐒[2, 2] * 𝐒[1, 2]'
    𝐙₃₂ = -𝐒[1, 1] * 𝐒[2, 1]' + 𝐒[2, 2] * 𝐒[1, 2]'
    𝐙₃₃ = 𝐒[1, 1] * 𝐒[2, 2]' + 𝐒[1, 2] * 𝐒[2, 1]'
    𝐙₃₄ = -1.0im * (𝐒[1, 1] * 𝐒[2, 2]' + 𝐒[2, 1] * 𝐒[1, 2]')
    𝐙₄₁ = 1.0im * (𝐒[2, 1] * 𝐒[1, 1]' + 𝐒[2, 2] * 𝐒[1, 2]')
    𝐙₄₂ = 1.0im * (𝐒[2, 1] * 𝐒[1, 1]' - 𝐒[2, 2] * 𝐒[1, 2]')
    𝐙₄₃ = -1.0im * (𝐒[2, 2] * 𝐒[1, 1]' - 𝐒[1, 2] * 𝐒[2, 1]')
    𝐙₄₄ = 𝐒[2, 2] * 𝐒[1, 1]' - 𝐒[1, 2] * 𝐒[2, 1]'

    𝐙 = @SMatrix [𝐙₁₁ 𝐙₁₂ 𝐙₁₃ 𝐙₁₄; 𝐙₂₁ 𝐙₂₂ 𝐙₂₃ 𝐙₂₄; 𝐙₃₁ 𝐙₃₂ 𝐙₃₃ 𝐙₃₄; 𝐙₄₁ 𝐙₄₂ 𝐙₄₃ 𝐙₄₄]
    return real.(𝐙)
end
