module EBCMPrecisionLossEstimators

using MLJ: machine, predict
using MLJXGBoostInterface
using TransitionMatrices: Spheroid, Cylinder, rmin, rmax

export estimate_integral_loss, estimate_inverse_loss, estimate_total_loss

for shape in [Spheroid, Cylinder]
    @eval function estimate_integral_loss(s::$shape, nₘₐₓ)
        mach = machine(joinpath(@__DIR__, "..", "fixtures", "$(string($shape)).Q.model"))
        predict(mach,
                [(rmax = rmax(s), rmin = rmin(s), nmax = nₘₐₓ, mr = real(s.m),
                  mi = imag(s.m))])[1]
    end

    @eval function estimate_inverse_loss(s::$shape, nₘₐₓ)
        mach = machine(joinpath(@__DIR__, "..", "fixtures", "$(string($shape)).T.model"))
        predict(mach,
                [(rmax = rmax(s), rmin = rmin(s), nmax = nₘₐₓ, mr = real(s.m),
                  mi = imag(s.m))])[1]
    end

    @eval function estimate_total_loss(s::$shape, nₘₐₓ)
        mach = machine(joinpath(@__DIR__, "..", "fixtures", "$(string($shape)).AT.model"))
        predict(mach,
                [(rmax = rmax(s), rmin = rmin(s), nmax = nₘₐₓ, mr = real(s.m),
                  mi = imag(s.m))])[1]
    end
end

end
