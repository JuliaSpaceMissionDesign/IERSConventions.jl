
 
export iers_polar_motion

""" 
    iers_polar_motion(m::IERSModel, xₚ::Number, yₚ::Number, tt_c::Number)

Compute the Polar Motion TIRF-to-ITRF rotation matrix, according to the IERS Conventions 
`m`, at time `tt_c` expressed in `TT` Julian centuries since `J2000`. The function requires 
`xp` and `yp`, the Celestial Intermediate Pole (CIP) coordinates with respect to the 
International Celestial Reference Frame (ITFR).

# TODO: expand description! (and check assumption for CPNd)

### References 
- IERS Technical Note No. [21](https://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn21.html)
- IERS Technical Note No. [32](https://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn32.html)
- IERS Technical Note No. [36](https://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn36.html) 
"""
function iers_polar_motion(m::IERSModel, xₚ::Number, yₚ::Number, tt_c::Number)
    sp = tio_locator(m, tt_c)
    return angle_to_dcm(sp, -xₚ, -yₚ, :ZYX)
end

iers_polar_motion(::CPND, xₚ::Number, yₚ::Number, ::Number) = DCM(1I)


"""
    tio_locator(m::IERSModel, tt_c::Number)

Compute the TIO locator `s'` at date, positioning the Terrestrial Intermediate Origin on 
the equator of the Celestial Intermediate Pole (CIP) at time `tt_c` expressed as `TT` Julian 
centuries since J2000. 

This function approximates the unpredictable motion of the TIO locator s' with its secular 
drift of ~0.47 μas/century. 

!!! note 
    A null value is returned for the [`iers1996`](@ref) and the [`CPNd`](@ref) models.

### References
- IERS Technical Note No. [36](https://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn36.html) 
- Lambert, S. and Bizouard C. (2002), Positioning the Terrestrial Ephemeris Origin in the 
  Terrestrial Reference Frame, [DOI: 10.1051/0004-6361:20021139](https://www.aanda.org/articles/aa/pdf/2002/40/aa2747.pdf)
"""
tio_locator(::IERSModel, tt_c::Number) = arcsec2rad(-47e-6tt_c) 

tio_locator(m::Union{IERS1996, CPND}, tt_c::Number) = 0 