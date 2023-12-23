
export iers_gmst, iers_gast

"""
    iers_gmst(m::IERSConventions, t::Number)

Compute the Greenwich Mean Sidereal Time (GMST), in radians, following the IERS Conventions 
`m` at time `t` expressed as `TT` Julian centuries since J2000. 

!!! note 
    The input time is automatically converted to UT1 for the computation of the Earth 
    Rotation Angle (ERA) or for the computation of the GMST of the 1996 conventions. Thus, 
    EOP data must be loaded before calling this function. 

### References 
- IERS Technical Note No. [21](https://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn21.html)
- IERS Technical Note No. [32](https://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn32.html)
- IERS Technical Note No. [36](https://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn36.html) 

### See also 
See also [`iers_gast`](@ref) and [`iers_era`](@ref).
"""
function iers_gmst(m::IERSConventions, t::Number)

    # Transform from TT centuries to UT1 days
    ut1 = t*Tempo.CENTURY2DAY + offset_tt2ut1(t*Tempo.CENTURY2SEC)/Tempo.DAY2SEC

    # Compute the Earth Rotation Angle 
    θ = orient_era(m, ut1)

    # Compute GMST 
    return iers_gmst(m, t, θ)

end


"""
    orient_gast(m::IERSConventions, t::Number, δΔψ::Number=0)

Compute the Greenwich Apparent Sidereal Time (GAST), in radians, following the IERS 
Conventions `m` at time `t` expressed as `TT` Julian centuries since J2000. An optional EOP 
correction for the nutation in longitude can be passed via `δΔψ` for the exact computation 
of the equation of the equinoxes. 

!!! note 
    The input time is automatically converted to UT1 for the computation of the Earth 
    Rotation Angle (ERA) or for the computation of the GMST of the 1996 conventions. Thus, 
    EOP data must be loaded before calling this function. 

### References 
- IERS Technical Note No. [21](https://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn21.html)
- IERS Technical Note No. [32](https://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn32.html)
- IERS Technical Note No. [36](https://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn36.html) 

### See also 
See also [`orient_gmst`](@ref), [`equation_equinoxes`](@ref) and [`orient_era`](@ref).
"""
function orient_gast(m::IERSConventions, t::Number, δΔψ::Number=0)
    return orient_gmst(m, t) + equation_equinoxes(m, t, δΔψ)
end


"""
    equation_equinoxes(m::IERSConventions, t::Number, δΔψ::Number = 0)

Compute the Equation of the Equinoxes, in radians, according to the IERS Conventions `m`, 
at time `t` expressed as `TT` Julian centuries since J2000. An optional EOP correction for 
the nutation in longitude can be passed via `δΔψ`.

!!! note 
    In this framework, the equation of the equinoxes is always defined according to the 
    following equation: `GAST = GMST + equation of the equinoxes`

### References 
- IERS Technical Note No. [21](https://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn21.html)
- IERS Technical Note No. [32](https://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn32.html)
- IERS Technical Note No. [36](https://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn36.html) 

### See also 
See also [`orient_obliquity`](@ref), [`orient_nutation_comp`](@ref) and [`eeq_complementary`](@ref).
"""
function equation_equinoxes(m::IERSConventions, t::Number, δΔψ::Number=0)

    # Retrive the mean obliquity and the nutation in longitude 
    ϵₐ = orient_obliquity(m, t)
    Δψ, _ = orient_nutation_comp(m, t) 

    return (Δψ + δΔψ)*cos(ϵₐ) + eeq_complementary(m, t)

end


"""
    eeq_complementary(m::IERSConventions, t::Number)

Compute the complementary terms of the equation of the equinoxes, in radians, associated 
to the IERS Conventions `m`, at time `t` expressed in `TT` Julian centuries since J2000.

!!! note 
    For the IERS 1996 Conventions, starting from 1997-02-27T00:00:00 UTC, two additional 
    terms are included to account for the Moon effect (see IAU 1997 C7 resolution). 
    The complementary terms of earlier dates are always null. 

### References 
- IAU Resolution C7, Recommendation 3 (1994)
- IERS Technical Note No. [21](https://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn21.html)
- IERS Technical Note No. [32](https://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn32.html)
- IERS Technical Note No. [36](https://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn36.html) 
"""
function eeq_complementary(m::IERSConventions, t::Number) 
    return ee_cpt(m, t, DelaunayArgs(m, t), PlanetaryArgs(m, t))
end


# 1996 CONVENTIONS
# ============================

# t = TT Julian centuries since J2000 
function iers_gmst(::IERS1996, t::Number)

    # Transform from TT centuries to UT1 centuries
    ut1 = t + offset_tt2ut1(t*Tempo.CENTURY2SEC)/Tempo.CENTURY2SEC

    # Coefficients for the IAU 1982 GMST-UT1 model. The first component has been adjusted 
    # of 12 hours because the input is defined with respect to noon of 01-01-2000
    A = -19089.45159
    B = 8640184.812866
    C = 0.093104
    D = -6.2e-6

    # Fractional part of UT1, in seconds
    f = Tempo.DAY2SEC*mod(ut1*Tempo.CENTURY2DAY, 1)

    # Compute GMST
    return mod2pi(2π/86400*(@evalpoly(ut1, A, B, C, D) + f))

end

function eeq_complementary(m::IERS1996, t::Number) 

    # Compute the longitude of the ascending node of the Moon
    Ω = delaunay_longitude_node(m, t)

    # Compute complementary terms 
    eeq = arcsec2rad(0.00264*sin(Ω) + 0.000063*sin(2Ω))

    # The IAU 1994 C7 resolution added two more terms to the equation of the equinoxes 
    # that are to be used from 1997-02-27T00:00:00 UTC   
    
    # The value of that date in TT Julian centuries since J2000 is: 
    t_ref = -8.9726337816e7/Tempo.CENTURY2SEC

    # Guarantees type stability
    return t < t_ref ? 0*eeq : eeq 

end


# 2003 CONVENTIONS
# ============================

function iers_gmst(::IERS2003, t::Number, θ::Number)

    # Evaluate the polynomial expression
    p = @evalpoly(
        t, 
        0.014506, 
     4612.15739966,   
        1.39667721, 
       -0.00009344,
        0.00001882
    )

    # Compute GMST 
    return mod2pi(θ + arcsec2rad(p))

end

include("constants/gast2000.jl")
build_series(:ee_cpt, :IERSConventions, [COEFFS_EEQ2000])


# 2010 CONVENTIONS
# ============================

# t = TT Julian centuries since J2000, θ = ERA 
function iers_gmst(::IERS2010, t::Number, θ::Number)

    # Evaluate interpolating series
    p = @evalpoly(
        t, 
        0.014506, 
     4612.156534,   
        1.3915817, 
       -0.00000044,
       -0.000029956, 
       -0.0000000368
    )

    # Compute GMST 
    return mod2pi(θ + arcsec2rad(p))
    
end