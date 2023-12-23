
export iers_obliquity

"""
    iers_obliquity(m::IERSConventions, t::Number)

Compute the mean obliquity of the ecliptic at epoch, in radians, at time `t` expressed in 
`TT` Julian centuries since `J2000`, according to the IERS convention `m`.

!!! note 
    The mean obliquity for the IERS 2003 conventions already accounts the adjustment to the 
    IAU 1976 precession model for the IAU 2000 precession-rate of the equator in obliquity .

### References 
- IERS Technical Note No. [21](https://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn21.html)
- IERS Technical Note No. [32](https://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn32.html)
- IERS Technical Note No. [36](https://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn36.html) 
"""
function iers_obliquity(::IERS1996, t::Number)
    return arcsec2rad(@evalpoly(t, 84381.448, -46.8150, -0.00059, 0.001813))
end

function iers_obliquity(::IERS2003, t::Number)
    δωₐ = -0.02524 # IAU 2003 Precession-rate adjustment
    return arcsec2rad(@evalpoly(t, 84381.448, -46.8150 + δωₐ, -0.00059, 0.001813))
end

function iers_obliquity(::IERS2010, t::Number)
    return arcsec2rad(        
        @evalpoly(
            t, 84381.406, -46.836769, -0.0001831, 0.00200340, -0.000000576, -0.0000000434
        )
    )
end