
"""
    fw_angles(m::IERS2010, t::Number) 

Compute the precession angles, γ, ϕ, ψ, ϵ in radians, according to the IAU 2006 
Fukushima-Williams 4-angle formulation at time `t` expressed in `TT` Julian centuries 
since `J2000`.

### References 
- Luzum, B. and Petit G. (2012), The IERS Conventions (2010), 
  [IERS Technical Note No. 36](https://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn36.html) 
- Wallace P. T. and Capitaine N. (2006), Precession-nutation procedures consistent with 
  IAU 2006 resolutions, [DOI: 10.1051/0004-6361:20065897](https://www.aanda.org/articles/aa/abs/2006/45/aa5897-06/aa5897-06.html) 
- [ERFA](https://github.com/liberfa/erfa/blob/master/src/pfw06.c) library
"""
function fw_angles(::IERS2010, t::Number)
    γ = arcsec2rad(
        @evalpoly(
            t, -0.052928, 10.556378, 0.4932044, -0.00031238, -0.000002788, 0.0000000260,
        )
    )

    ϕ = arcsec2rad(
        @evalpoly(
            t, 84381.412819, -46.811016, 0.0511268, 0.00053289, -0.000000440, -0.0000000176,
        )
    )

    ψ = arcsec2rad(
        @evalpoly(
            t, -0.041775, 5038.481484, 1.5584175, -0.00018522, -0.000026452, -0.0000000148,
        )
    )

    ϵ = iers_obliquity(iers2010a, t)
    
    return γ, ϕ, ψ, ϵ
end


"""
    fw2xy(ϵ::Number, ψ::Number, γ::Number, φ::Number)

Compute the CIP X and Y coordinates from Fukushima-Williams bias-precession-nutation 
angles, in radians.

### References
- Wallace P. T. and Capitaine N. (2006), Precession-nutation procedures consistent with 
  IAU 2006 resolutions, [DOI: 10.1051/0004-6361:20065897](https://www.aanda.org/articles/aa/abs/2006/45/aa5897-06/aa5897-06.html) 
"""
function fw2xy(γ::Number, φ::Number, ψ::Number, ϵ::Number)

    sϵ, cϵ = sincos(ϵ)
    sψ, cψ = sincos(ψ)
    sγ, cγ = sincos(γ)
    sϕ, cϕ = sincos(φ)

    a = (sϵ * cψ * cϕ - cϵ * sϕ)

    x = sϵ * sψ * cγ - a * sγ
    y = sϵ * sψ * sγ + a * cγ

    return x, y
end


"""
    fw_matrix(γ, ϕ, ψ, ε)

Form a rotation matrix given the Fukushima-Williams angles, expressed in radians.

The present function can construct three different matrices depending on which angles are 
supplied as arguments: 

- **NPB**: To obtain the Nutation-Precession-Bias (NPB) matrix, generate the four standard FW 
    precession angles (̄γ, ̄ϕ, ̄ψ, ϵₐ) then generate the nutation components Δψ and Δϵ and add them 
    to ̄ψ, ϵₐ. Finally, call the present functions using those four angles as arguments. 

- **PB**: To obtain the Precession-Bias matrix, generate the four standard FW precession 
    angles and call the present function. 

- **B**: To obtain the frame bias matrix, generate the four standard FW precession angles at 
    date J2000.0 and call this function.

The remaining nutation-only and precession only matrices can be obtained by combining these 
three appropriately. 

### References
- Wallace P. T. and Capitaine N. (2006), Precession-nutation procedures consistent with 
  IAU 2006 resolutions, [DOI: 10.1051/0004-6361:20065897](https://www.aanda.org/articles/aa/abs/2006/45/aa5897-06/aa5897-06.html) 
- [ERFA](https://github.com/liberfa/erfa/blob/master/src/fw2m.c) library
"""
function fw_matrix(γ, ϕ, ψ, ϵ)
    return angle_to_dcm(-ϵ, :X) * angle_to_dcm(γ, ϕ, -ψ, :ZXZ)
end