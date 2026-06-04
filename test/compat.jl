@testitem "Numeric compatibility shims" begin
    using TransitionMatrices
    using Arblib
    using ForwardDiff

    @test precision(ComplexF32) == 24
    @test precision(ComplexF64) == 53
    @test precision(Complex{Double64}) == 106
    @test precision(ComplexF128) == 113
    @test precision(Complex{Arb}) == precision(Arb)

    Tag = ForwardDiff.Tag{typeof(identity), Float64}
    DualT = ForwardDiff.Dual{Tag, Float64, 1}
    @test precision(Complex{DualT}) == precision(Float64)

    arb = Arb(1.25)
    @test convert(Float128, arb) == Float128(1.25)
    @test TransitionMatrices.Quadmath.Float128(arb) == Float128(1.25)
    @test Float32(arb) == Float32(1.25)
    @test Int64(Arb(2.0)) == 2
    @test round(Arb(1.2), RoundUp) == ceil(BigFloat(Arb(1.2)))
    @test complex(Arb) === Acb

    dual = ForwardDiff.Dual(1.5, 2.0)
    @test Float16(dual) == Float16(1.5)
    @test Float32(dual) == Float32(1.5)
    @test BigFloat(dual) == BigFloat(1.5)

    target = Arb()
    Arblib.set!(target, dual)
    @test Float64(target) == 1.5
    Arblib.set!(target, Float128(1.75))
    @test Float64(target) == 1.75
    Arblib.set!(target, Double64(2.25))
    @test Float64(target) == 2.25

    acb = Acb(1, 2)
    @test abs2(acb) == Arb(5)
    converted_arb = convert(Complex{Arb}, acb)
    @test real(converted_arb) == Arb(1)
    @test imag(converted_arb) == Arb(2)

    converted_dual = convert(Complex{DualT}, acb)
    @test ForwardDiff.value(real(converted_dual)) == 1.0
    @test ForwardDiff.value(imag(converted_dual)) == 2.0

    A = Arb[1 2; 3 5]
    Ainv = inv(A)
    ArbI = A * Ainv
    @test Float64(abs(real(ArbI[1, 1] - 1))) < 1e-30
    @test Float64(abs(real(ArbI[1, 2]))) < 1e-30
    @test Float64(abs(real(ArbI[2, 1]))) < 1e-30
    @test Float64(abs(real(ArbI[2, 2] - 1))) < 1e-30

    C = Acb[1 + 0im 2 + 0im; 3 + 0im 5 + 0im]
    Cinv = inv(C)
    AcbI = C * Cinv
    @test Float64(abs(real(AcbI[1, 1] - 1))) < 1e-30
    @test Float64(abs(real(AcbI[1, 2]))) < 1e-30
    @test Float64(abs(real(AcbI[2, 1]))) < 1e-30
    @test Float64(abs(real(AcbI[2, 2] - 1))) < 1e-30
end
