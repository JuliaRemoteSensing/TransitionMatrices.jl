module RadiativeTransferShared

using NCDatasets: defVar

export R_DRY_AIR, RAIN_N0, WATER_DENSITY, MASS_COEFF, air_density, write_var

const R_DRY_AIR = 287.15
const RAIN_N0 = 8.0e3
const WATER_DENSITY = 997.0
const MASS_COEFF = pi * WATER_DENSITY / 6 * 1e-9

function air_density(qvapor, pb, p, theta_perturbation)
    pressure = pb + p
    tempk = (theta_perturbation + 300.0) * (pressure / 100000.0)^(2.0 / 7.0)
    0.622 * pressure / (R_DRY_AIR * tempk * (qvapor + 0.622))
end

function write_var(ds, name, data, dims; units = "", description = "")
    v = defVar(ds, name, data, dims)
    isempty(units) || (v.attrib["units"] = units)
    isempty(description) || (v.attrib["description"] = description)
    return v
end

end
