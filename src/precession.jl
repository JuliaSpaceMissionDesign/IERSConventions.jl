
export iers_precession

""" 
    iers_precession(m::IERSModel, t::Number)

Return the precession matrix that rotates a vector from MEME2000 axes to Mean of Date (MOD) 
axes, at time `t` expressed in `TT` Julian centuries since `J2000`, according to the IERS 
convention `m`. 

!!! note    
    This matrix rotates vectors from the Mean Equator and Mean Equinox of J2000 (MEME2000) 
    to Mean-of-Date (MOD) axes. The frame bias between the GCRF and MEME2000 is excluded 
    from the returned matrix and must eventually be included with a separate rotation.

!!! note 
    The IAU Working Group on Precession and the Ecliptic (Hilton, 2006) has decided to leave 
    the choice of the parameterization for the precession angles to the user. In this 
    function for the IERS 1996 conventions, the precession matrix is computed using the 
    traditional parameterization of Newcomb and Liekse (`zₐ`, `θₐ`, `ζₐ`), whereas it adopts 
    the 4-angles formulation (`ϵ₀`, `ψₐ`, `ωₐ`, `χₐ`) recommended by (Capitaine et al., 2003a)
    for all the remaining models.

### References 
- Lieske J. H. et al, (1977), Expression for the Precession Quantities Based upon the IAU 
    (1976) System of Astronomical Constants.
- Hilton J. L. et al., (2006), Report of the International Astronomical Union Division I Working 
    Group on Precession and the Ecliptic.
- Capitaine N. et al., (2003a), Expressions for IAU 2000 precession quantities. 
- IERS Technical Note No. [21](https://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn21.html)
- IERS Technical Note No. [32](https://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn32.html)
- IERS Technical Note No. [36](https://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn36.html) 

### See also 
See also [`precession_angles_rot3`](@ref) and [`precession_angles_rot4`](@ref)
"""
function iers_precession(m::IERSModel, t::Number)

    # Compute the precession angles 
    ϵ₀, ψₐ, ωₐ, χₐ = precession_angles_rot4(m, t)

    # Form the precession matrix
    return angle_to_dcm(χₐ, :Z)*angle_to_dcm(ϵ₀, -ψₐ, -ωₐ, :XZX)

end


""" 
    precession_angles_rot3(m::IERSModel, t::Number)

Return the equatorial precession angles zₐ, θₐ, ζₐ, in radians, at time `t` expressed in 
`TT` Julian centuries since `J2000` for the 3-rotations precession series initially used 
by Newcomb and Lieske. 

!!! note 
    The expressions for these angles compatible with the IAU 2000A precession and nutation 
    have been developed in order to match the 4-rotation series to sub-microarcsecond 
    accuracy over 4 centuries.

!!! note 
    `TT` is used as time argument instead of TDB. The largest term in the TDB-TT difference 
    causes a periodic error which is significantly under the microarcsecond accuracy.

### References 
- Lieske J. H. et al, (1977), Expression for the Precession Quantities Based upon the IAU 
    (1976) System of Astronomical Constants.
- Capitained N. 2003b et al., Expressions for IAU 2000 precession quantities.
- Hilton J. L. et al., (2006), Report of the International Astronomical Union Division I Working 
Group on Precession and the Ecliptic.
- IERS Technical Note No. [21](https://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn21.html)
- IERS Technical Note No. [32](https://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn32.html)
- IERS Technical Note No. [36](https://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn36.html) 

### See also 
See also [`iers_precession`](@ref) and [`precession_angles_rot4`](@ref)
"""
precession_angles_rot3


"""

    precession_angles_rot4(m::IERSModel, t::Number)

Return the precession angles ϵ₀, ψₐ, ωₐ, χₐ, in radians, at time `t` expressed in `TT` Julian 
centuries since `J2000` required for the canonical 4-rotations precession series.

!!! note 
    `TT` is used as time argument instead of TDB. The largest term in the TDB-TT difference 
    causes a periodic error which is significantly under the microarcsecond accuracy.

### References 
- IERS Technical Note No. [21](https://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn21.html)
- IERS Technical Note No. [32](https://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn32.html)
- IERS Technical Note No. [36](https://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn36.html) 

### See also 
See also [`iers_precession`](@ref) and [`precession_angles_rot4`](@ref)
"""
precession_angles_rot4


# 1996 CONVENTIONS
# ============================

function iers_precession(m::IERS1996, t::Number)

    # Compute the precession angles 
    zₐ, θₐ, ζₐ = precession_angles_rot3(m, t)

    # Form the precession matrix
    return angle_to_dcm(-ζₐ, θₐ, -zₐ, :ZYZ)

end

function precession_angles_rot3(::IERS1996, t::Number)

    ζₐ = arcsec2rad(@evalpoly(t, 0, 2306.2181,  0.30188,  0.017998))
    θₐ = arcsec2rad(@evalpoly(t, 0, 2004.3109, -0.42665, -0.041833))
    zₐ = arcsec2rad(@evalpoly(t, 0, 2306.2181,  1.09468,  0.018203))

    return zₐ, θₐ, ζₐ

end

function precession_angles_rot4(m::IERS1996, t::Number)

    ϵ₀ = iers_obliquity(m, 0)
    ωₐ = ϵ₀ + arcsec2rad(@evalpoly(t, 0, 0, 0.05127, -0.007726))

    ψₐ = arcsec2rad(@evalpoly(t, 0, 5038.7784, -1.07259, -0.001147))
    χₐ = arcsec2rad(@evalpoly(t, 0, 10.5526, -2.38064, -0.001125))

    return ϵ₀, ψₐ, ωₐ, χₐ

end


# 2003 CONVENTIONS
# ============================

function precession_angles_rot3(::IERS2003, t::Number)

    ζₐ = arcsec2rad(
        @evalpoly(t, 2.5976176, 2306.0809506,  0.3019015, 0.0179663, -0.0000327, -0.0000002)
    )

    θₐ = arcsec2rad(
        @evalpoly(t, 0, 2004.1917476, -0.4269353, -0.0418251, -0.0000601, -0.0000001)
    )

    zₐ = arcsec2rad(
        @evalpoly(t, -2.5976176, 2306.0803226, 1.0947790,  0.0182273, 0.0000470, -0.0000003)
    )

    return zₐ, θₐ, ζₐ

end

function precession_angles_rot4(m::IERS2003, t::Number)

    ϵ₀ = iers_obliquity(m, 0)
    ωₐ = ϵ₀ + arcsec2rad(@evalpoly(t, 0, -0.02524, 0.05127, -0.007726))

    ψₐ = arcsec2rad(@evalpoly(t, 0, 5038.47875, -1.07259, -0.001147))
    χₐ = arcsec2rad(@evalpoly(t, 0, 10.5526, -2.38064, -0.001125))

    return ϵ₀, ψₐ, ωₐ, χₐ

end


# 2010 CONVENTIONS
# ============================

function precession_angles_rot3(::IERS2010, t::Number)

    # Values have been retrieved from Table 1 of Hilton (2006)
    ζₐ = arcsec2rad(
        @evalpoly(t, 2.650545, 2306.083227, 0.2988499, 0.01801828, -5.971e-6, -3.173e-7)
    )

    θₐ = arcsec2rad(
        @evalpoly(t, 0, 2004.191903, -0.4294934, -0.04182264, -7.089e-6, -1.274e-7)
    )

    zₐ = arcsec2rad(
        @evalpoly(t, -2.650545, 2306.077181, 1.0927348, 0.01826837, -0.000028596, -2.904e-7)
    )

    return zₐ, θₐ, ζₐ

end

function precession_angles_rot4(m::IERS2010, t::Number)
    
    # This is the so-called P03 model!
    ϵ₀ = iers_obliquity(m, 0)

    ωₐ = ϵ₀ + arcsec2rad(
        @evalpoly(t, 0, -0.025754, 0.0512623, -0.00772503, -0.000000467, 0.0000003337)
    )

    ψₐ = arcsec2rad(
        @evalpoly(t, 0, 5038.481507, -1.0790069, -0.00114045, 0.000132851, -0.0000000951)
    )

    χₐ = arcsec2rad(
        @evalpoly(t, 0, 10.556403, -2.3814292, -0.00121197, 0.000170663, -0.0000000560)
    )

    return ϵ₀, ψₐ, ωₐ, χₐ

end


