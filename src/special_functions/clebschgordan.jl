function clebschgordan(j₁, m₁, j₂, m₂, j₃, m₃ = m₁ + m₂)
    if !(j₃ <= j₁ + j₂ && j₁ <= j₂ + j₃ && j₂ <= j₃ + j₁)
        return 0.0
    end

    return (-1.0)^(j₁ - j₂ + m₃) * √(2j₃ + 1) * wig3jj(2j₁, 2j₂, 2j₃, 2m₁, 2m₂, -2m₃)
end

@testitem "Clebsch-Gordan coefficients should be correct" begin
    using Wigxjpf: wig_table_init, wig_table_free, wig_temp_init, wig_temp_free
    using WignerSymbols: WignerSymbols

    wig_table_init(10, 3)
    wig_temp_init(10)

    for j₁ in 0:5, m₁ in (-j₁):j₁, j₂ in 0:5, m₂ in (-j₂):j₂, j₃ in 0:5, m₃ in (-j₃):j₃
        @test TransitionMatrices.clebschgordan(j₁, m₁, j₂, m₂, j₃, m₃) ≈
              WignerSymbols.clebschgordan(j₁, m₁, j₂, m₂, j₃, m₃)
    end

    wig_temp_free()
    wig_table_free()
end
