export  iers_rot3_gcrf_to_mod, 
        iers_rot3_gcrf_to_tod, 
        iers_rot3_gcrf_to_gtod, 
        iers_rot3_gcrf_to_pef, 
        iers_rot3_itrf_to_mod, 
        iers_rot3_itrf_to_tod, 
        iers_rot3_itrf_to_gtod, 
        iers_rot3_itrf_to_pef, 
        iers_rot3_gcrf_to_cirf, 
        iers_rot3_gcrf_to_itrf,
        iers_rot3_gcrf_to_tirf, 
        iers_rot3_itrf_to_cirf, 
        iers_rot3_itrf_to_tirf,

        iers_rot6_gcrf_to_itrf,
        iers_rot9_gcrf_to_itrf,
        iers_rot12_gcrf_to_itrf

# Equinox-based transformations from GCRF
# ==============================================

"""
    iers_rot3_gcrf_to_mod(tt_s::Number, m::IERSModel=iers2010b)

Compute the rotation matrix from the Geocentric Celestial Reference Frame (GCRF) to 
the Mean-of-Date (MOD) at time `tt_s`, expressed in TT seconds since `J2000`.

!!! note 
    The Mean-of-Date axes are obtained by applying the frame bias and precession matrix. 
    For this reason, if the [`iers1996`](@ref) conventions are used, the rotation is 
    actually computed starting from the MEME2000 rather than the GCRF.  

### See also 
See also [`iers_pb`](@ref) and [`iers_rot3_itrf_to_mod`](@ref).

"""
function iers_rot3_gcrf_to_mod(tt_s::Number, m::IERSModel=iers2010b)
    return iers_pb(m, tt_s/Tempo.CENTURY2SEC)
end


"""
    iers_rot3_gcrf_to_tod(tt_s::Number, m::IERSModel=iers2010b)

Compute the rotation matrix from the Geocentric Celestial Reference Frame (GCRF) to 
the True-of-Date (TOD) at time `tt_s`, expressed in TT seconds since `J2000`.

!!! note 
    The True-of-Date axes are obtained by applying the frame bias, precession and 
    nutation matrix. For this reason, if the [`iers1996`](@ref) conventions are used, the 
    rotation is actually computed starting from the MEME2000 rather than the GCRF.  

!!! note 
    The EOP nutation corrections are only used for the [`iers2003a`](@ref) and 
    [`iers2010a`](@ref) models. 

## See also 
See also [`iers_npb`](@ref) and [`iers_rot3_itrf_to_tod`](@ref).
"""
function iers_rot3_gcrf_to_tod(tt_s::Number, m::IERSModel=iers2010b)

    ttc = tt_s/Tempo.CENTURY2SEC

    # Retrieve the EOP corrections to the nutation in longitude and obliquity
    δΔψ = eop_δΔψ(m, ttc)
    δΔϵ = eop_δΔϵ(m, ttc)

    # Compute Nutation-Precession-Bias matrix
    return iers_npb(m, ttc, δΔψ, δΔϵ)

end


"""
    iers_rot3_gcrf_to_gtod(tt_s::Number, m::IERSModel=iers2010b)

Compute the rotation matrix from the Geocentric Celestial Reference Frame (GCRF) to 
the Greenwich True-of-Date (GTOD) at time `tt_s`, expressed in TT seconds since `J2000`.

!!! note 
    If the [`iers1996`](@ref) conventions are used, the rotation is actually computed 
    starting from the MEME2000 rather than the GCRF.  

!!! note 
    The EOP nutation corrections are only used for the [`iers2003a`](@ref) and 
    [`iers2010a`](@ref) models. 

## See also 
See also [`iers_rot3_itrf_to_gtod`](@ref).
"""
function iers_rot3_gcrf_to_gtod(tt_s::Number, m::IERSModel=iers2010b)

    ttc = tt_s/Tempo.CENTURY2SEC

    # Retrieve the EOP corrections to the nutation in longitude and obliquity
    δΔψ = eop_δΔψ(m, ttc)
    δΔϵ = eop_δΔϵ(m, ttc)

    # Retrieve precession-bias matrix
    PB = iers_pb(m, ttc)

    # Retrieve MOD-to-GTOD rotation 
    RN = mod_to_gtod3(ttc, m, δΔψ, δΔϵ)
    return RN*PB

end


"""
    iers_rot3_gcrf_to_pef(tt_s::Number, m::IERSModel=iers2010b)

Compute the rotation matrix from the Geocentric Celestial Reference Frame (GCRF) to 
the Pseudo-Earth Fixed (PEF) at time `tt_s`, expressed in TT seconds since `J2000`.

!!! note 
    If the [`iers1996`](@ref) conventions are used, the rotation is actually computed 
    starting from the MEME2000 rather than the GCRF.  

!!! note 
    The EOP nutation corrections are only used for the [`iers2003a`](@ref) and 
    [`iers2010a`](@ref) models. 

!!! note 
    For the [`iers1996`](@ref) and [`CPNd`](@ref) models, there are no differences 
    between the PEF and GTOD axes. # TODO: specify the magnitude of the rotation!

## See also 
See also [`iers_rot3_gcrf_to_gtod`](@ref) and [`iers_rot3_itrf_to_pef`](@ref).
"""
function iers_rot3_gcrf_to_pef(tt_s::Number, m::IERSModel=iers2010b)

    # Compute the GCRF to GTOD rotation matrix
    RNPB = iers_rot3_gcrf_to_gtod(tt_s, m)

    # Compute the TIO locator matrix (from GTOD to PEF)
    sp = tio_locator(m, tt_s/Tempo.CENTURY2SEC)

    # Assemble the global rotation
    return angle_to_dcm(sp, :Z)*RNPB

end

# TODO: add TEME rotations


# Equinox-based transformations from ITRF
# ==============================================

"""
    iers_rot3_itrf_to_pef(tt_s::Number, m::IERSModel=iers2010b)

Compute the rotation matrix from the International Terrestrial Reference Frame (ITRF) to 
the Pseudo-Earth Fixed (PEF) at time `tt_s`, expressed in TT seconds since `J2000`.

!!! note 
    The [`CPNd`](@ref) model returns an identity rotation because it neglects the 
    effects of polar motion. # TODO: specify magnitude!

## See also 
See also [`iers_rot3_gcrf_to_pef`](@ref).
"""
function iers_rot3_itrf_to_pef(tt_s::Number, m::IERSModel=iers2010b)

    ttc = tt_s/Tempo.CENTURY2SEC

    # Retrieve pole coordinates 
    xₚ = eop_xp(m, ttc)
    yₚ = eop_yp(m, ttc)

    # Compute the polar motion rotation without the TIO contribution 
    # Note: this equation is wrongly reported on Vallado Fundamentals of Astrodyn. 4th Ed. 
    return angle_to_dcm(yₚ, xₚ, :XY)

end

iers_rot3_itrf_to_pef(::Number, ::CPND) = DCM(1I)


"""
    iers_rot3_itrf_to_gtod(tt_s::Number, m::IERSModel=iers2010b)

Compute the rotation matrix from the International Terrestrial Reference Frame (ITRF) to 
the Greenwich True-of-Date (GTOD) at time `tt_s`, expressed in TT seconds since `J2000`.

!!! note 
    The [`CPNd`](@ref) model returns an identity rotation because it neglects the 
    effects of polar motion. # TODO: specify magnitude!

!!! note 
    For the [`iers1996`](@ref) model, the GTOD and PEF axes are equal because the IERS 
    1996 conventions do not account for the TIO locator effects.

## See also 
See also [`iers_rot3_itrf_to_pef`](@ref) and [`iers_rot3_gcrf_to_gtod`](@ref).
"""
function iers_rot3_itrf_to_gtod(tt_s::Number, m::IERSModel=iers2010b)

    # Compute ITRF to PEF rotation
    W  = iers_rot3_itrf_to_pef(tt_s, m)

    # Compute TIO locator rotation 
    sp = tio_locator(m, tt_s/Tempo.CENTURY2SEC)
    return angle_to_dcm(-sp, :Z)*W

end


"""
    iers_rot3_itrf_to_tod(tt_s::Number, m::IERSModel=iers2010b)

Compute the rotation matrix from the International Terrestrial Reference Frame (ITRF) to 
the True-of-Date (TOD) at time `tt_s`, expressed in TT seconds since `J2000`.

!!! note 
    The EOP corrections in longitude are only used for the [`iers2003a`](@ref) and 
    [`iers2010a`](@ref) models. 

## See also 
See also [`iers_rot3_itrf_to_gtod`](@ref) and [`iers_rot3_gcrf_to_tod`](@ref).
"""
function iers_rot3_itrf_to_tod(tt_s::Number, m::IERSModel=iers2010b)

    ttc = tt_s/Tempo.CENTURY2SEC

    # Compute EOP nutation corrections 
    δΔψ = eop_δΔψ(m, ttc)

    # Compute ITRF to GTOD rotation 
    W = iers_rot3_itrf_to_gtod(tt_s, m)

    # Compute GTOD to TOD rotation matrix 
    R = angle_to_dcm(-iers_gast(m, ttc, δΔψ), :Z)
    return R*W

end


"""
    iers_rot3_itrf_to_mod(tt_s::Number, m::IERSModel=iers2010b)

Compute the rotation matrix from the International Terrestrial Reference Frame (ITRF) to 
the Mean-of-Date (MOD) at time `tt_s`, expressed in TT seconds since `J2000`.

!!! note 
    The EOP nutation corrections are only used for the [`iers2003a`](@ref) and 
    [`iers2010a`](@ref) models. 

## See also 
See also [`iers_rot3_itrf_to_tod`](@ref) and [`iers_rot3_gcrf_to_mod`](@ref).
"""
function iers_rot3_itrf_to_mod(tt_s::Number, m::IERSModel=iers2010b)

    ttc = tt_s/Tempo.CENTURY2SEC

    # Compute ITRF to GTOD rotation 
    W = iers_rot3_itrf_to_gtod(tt_s, m)

    # Retrieve the EOP corrections to the nutation in longitude and obliquity
    δΔψ = eop_δΔψ(m, ttc)
    δΔϵ = eop_δΔϵ(m, ttc)

    # Compute GTOD-to-MOD rotation
    RNt = mod_to_gtod3(ttc, m, δΔψ, δΔϵ)' 
    return RNt*W 

end


# CIO-based transformations from GCRF 
# ==============================================

"""
    iers_rot3_gcrf_to_cirf(tt_s::Number, m::IERSModel=iers2010b)

Compute the rotation matrix from the Geocentric Celestial Reference Frame (GCRF) to 
the Celestial Intermediate Reference Frame (CIRF) at time `tt_s`, expressed in TT seconds 
since `J2000`, following the IERS Conventions `m`. 

!!! note 
    EOP corrections to the CIP coordinates (δX, δY) are only added in the [`iers2003a`](@ref)
    and [`iers2010a`](@ref) models.

### References 
- IERS Technical Note No. [36](https://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn36.html) 

### See also 
See also [`iers_rot3_gcrf_to_tirf`](@ref) and [`iers_rot3_gcrf_to_itrf`](@ref).
"""
function iers_rot3_gcrf_to_cirf(tt_s::Number, m::IERSModel=iers2010b)    
    ttc = tt_s/Tempo.CENTURY2SEC
    return iers_cip_motion(m, ttc, eop_δX(m, ttc), eop_δY(m, ttc))
end

"""
    iers_rot3_gcrf_to_tirf(tt_s::Number, m::IERSModel=iers2010b)

Compute the rotation matrix from the Geocentric Celestial Reference Frame (GCRF) to 
the Terrestrial Intermediate Reference Frame (TIRF) at time `tt_s`, expressed in TT seconds 
since `J2000`, following the IERS Conventions `m`. 

!!! note 
    EOP corrections to the CIP coordinates (δX, δY) are only added in the [`iers2003a`](@ref)
    and [`iers2010a`](@ref) models.

### References 
- IERS Technical Note No. [36](https://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn36.html) 

### See also 
See also [`iers_rot3_gcrf_to_cirf`](@ref) and [`iers_rot3_gcrf_to_itrf`](@ref).
"""
function iers_rot3_gcrf_to_tirf(tt_s::Number, m::IERSModel=iers2010b)

    # Convert TT seconds to UT1 days since J2000
    ut1_d = (tt_s + offset_tt2ut1(tt_s))/Tempo.DAY2SEC

    # Compute GCRF to CIRF rotation (CIP motion)
    Q = iers_rot3_gcrf_to_cirf(tt_s, m)

    # Compute the ERA rotation matrixes
    R = iers_era_rotm(m, ut1_d)

    return R*Q

end


"""
    iers_rot3_gcrf_to_itrf(tt_s::Number, m::IERSModel=iers2010b)

Compute the rotation matrix from the Geocentric Celestial Reference Frame (GCRF) to 
the International Terrestrial Reference Frame (ITRF) at time `tt_s`, expressed in TT seconds 
since `J2000`, following the IERS Conventions `m`. 

!!! note 
    EOP corrections to the CIP coordinates (δX, δY) are only added in the [`iers2003a`](@ref)
    and [`iers2010a`](@ref) models.

!!! note 
    Polar motion is neglected in the [`CPNd`](@ref) model. 

### References 
- IERS Technical Note No. [36](https://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn36.html) 

### See also 
See also [`iers_rot3_gcrf_to_cirf`](@ref) and [`iers_rot3_gcrf_to_tirf`](@ref).
"""
function iers_rot3_gcrf_to_itrf(tt_s::Number, m::IERSModel=iers2010b)

    # Compute the GCRF to TIRF rotation matrix
    RQ = iers_rot3_gcrf_to_tirf(tt_s, m)

    ttc = tt_s/Tempo.CENTURY2SEC

    # Retrieve pole coordinates 
    xₚ = eop_xp(m, ttc)
    yₚ = eop_yp(m, ttc)

    # Compute the polar motion rotation
    W = iers_polar_motion(m, xₚ, yₚ, ttc)

    return W*RQ

end

"""
    iers_rot6_gcrf_to_itrf(tt_s::Number, m::IERSModel=iers2010b)

Compute the rotation matrix and its derivative from the Geocentric Celestial Reference Frame 
(GCRF) to the International Terrestrial Reference Frame (ITRF) at time `tt_s`, expressed 
in TT seconds since `J2000`, following the IERS Conventions `m`. 

!!! note 
    EOP corrections to the CIP coordinates (δX, δY) are only added in the [`iers2003a`](@ref)
    and [`iers2010a`](@ref) models.

!!! note 
    Polar motion is neglected in the [`CPNd`](@ref) model. 

!!! note 
    The time derivative of the nutation and precession effects is neglected.

### References 
- IERS Technical Note No. [36](https://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn36.html) 
"""
function iers_rot6_gcrf_to_itrf(tt_s::Number, m::IERSModel=iers2010b)

    RQ = iers_rot3_gcrf_to_tirf(tt_s, m)

    ttc = tt_s/Tempo.CENTURY2SEC

    # Retrieve pole coordinates 
    xₚ = eop_xp(m, ttc)
    yₚ = eop_yp(m, ttc)
    LOD = eop_lod(m, ttc)

    ωe = SVector(0.0, 0.0, iers_earth_rot_rate(LOD))
    Ω = -skew(ωe)

    # Compute the polar motion rotation
    W = iers_polar_motion(m, xₚ, yₚ, ttc)

    return W*RQ, W*Ω*RQ

end

"""
    iers_rot9_gcrf_to_itrf(tt_s::Number, m::IERSModel=iers2010b)

Compute the rotation matrix, its first and second derivative from the Geocentric Celestial 
Reference Frame (GCRF) to the International Terrestrial Reference Frame (ITRF) at time 
`tt_s`, expressed in TT seconds since `J2000`, following the IERS Conventions `m`. 

!!! note 
    EOP corrections to the CIP coordinates (δX, δY) are only added in the [`iers2003a`](@ref)
    and [`iers2010a`](@ref) models.

!!! note 
    Polar motion is neglected in the [`CPNd`](@ref) model. 

!!! note 
    The time derivative of the nutation and precession effects is neglected.

### References 
- IERS Technical Note No. [36](https://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn36.html) 
"""
function iers_rot9_gcrf_to_itrf(tt_s::Number, m::IERSModel=iers2010b)

    RQ = iers_rot3_gcrf_to_tirf(tt_s, m)

    ttc = tt_s/Tempo.CENTURY2SEC

    # Retrieve pole coordinates 
    xₚ = eop_xp(m, ttc)
    yₚ = eop_yp(m, ttc)
    LOD = eop_lod(m, ttc)

    ωe = SVector(0.0, 0.0, iers_earth_rot_rate(LOD))
    Ω = -skew(ωe)

    # Compute the polar motion rotation
    W = iers_polar_motion(m, xₚ, yₚ, ttc)

    ΩRW = Ω*RQ

    return W*RQ, W*ΩRW, W*Ω*ΩRW
end

"""
    iers_rot12_gcrf_to_itrf(tt_s::Number, m::IERSModel=iers2010b)

Compute the rotation matrix, its first, second and third derivative from the Geocentric 
Celestial Reference Frame (GCRF) to the International Terrestrial Reference Frame (ITRF) 
at time `tt_s`, expressed in TT seconds since `J2000`, following the IERS Conventions `m`. 

!!! note 
    EOP corrections to the CIP coordinates (δX, δY) are only added in the [`iers2003a`](@ref)
    and [`iers2010a`](@ref) models.

!!! note 
    Polar motion is neglected in the [`CPNd`](@ref) model. 

!!! note 
    The time derivative of the nutation and precession effects is neglected.

### References 
- IERS Technical Note No. [36](https://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn36.html) 
"""
function iers_rot12_gcrf_to_itrf(tt_s::Number, m::IERSModel=iers2010b)

    RQ = iers_rot3_gcrf_to_tirf(tt_s, m)

    ttc = tt_s/Tempo.CENTURY2SEC

    # Retrieve pole coordinates 
    xₚ = eop_xp(m, ttc)
    yₚ = eop_yp(m, ttc)
    LOD = eop_lod(m, ttc)

    ωe = SVector(0.0, 0.0, iers_earth_rot_rate(LOD))
    Ω = -skew(ωe)

    # Compute the polar motion rotation
    W = iers_polar_motion(m, xₚ, yₚ, ttc)

    ΩRW = Ω*RQ
    Ω²RW = Ω*ΩRW

    return W*RQ, W*ΩRW, W*Ω²RW, W*Ω*Ω²RW
end


# CIO-based transformations from ITRF 
# ==============================================

"""
    iers_rot3_itrf_to_tirf(tt_s::Number, m::IERSModel=iers2010b)

Compute the rotation matrix from the International Terrestrial Reference Frame (ITRF) to 
the Terrestrial Intermediate Reference Frame (TIRF) at time `tt_s`, expressed in TT seconds 
since `J2000`, following the IERS Conventions `m`. 

!!! note 
    Polar motion is neglected in the [`CPNd`](@ref) model. 

### References 
- IERS Technical Note No. [36](https://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn36.html) 

### See also 
See also [`iers_rot3_gcrf_to_tirf`](@ref) and [`iers_rot3_itrf_to_cirf`](@ref).
"""
function iers_rot3_itrf_to_tirf(tt_s::Number, m::IERSModel=iers2010b)

    ttc = tt_s/Tempo.CENTURY2SEC

    # Retrieve pole coordinates 
    xₚ = eop_xp(m, ttc)
    yₚ = eop_yp(m, ttc)

    # Compute the polar motion rotation
    return iers_polar_motion(m, xₚ, yₚ, ttc)'

end

iers_rot3_itrf_to_tirf(::Number, ::CPND) = DCM(1I)


"""
    iers_rot3_itrf_to_cirf(tt_s::Number, m::IERSModel=iers2010b)

Compute the rotation matrix from the International Terrestrial Reference Frame (ITRF) to 
the Celestial Intermediate Reference Frame (CIRF) at time `tt_s`, expressed in TT seconds 
since `J2000`, following the IERS Conventions `m`. 

!!! note 
    Polar motion is neglected in the [`CPNd`](@ref) model. 

### References 
- IERS Technical Note No. [36](https://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn36.html) 

### See also 
See also [`iers_rot3_gcrf_to_cirf`](@ref) and [`iers_rot3_itrf_to_tirf`](@ref).
"""
function iers_rot3_itrf_to_cirf(tt_s::Number, m::IERSModel=iers2010b)

    # Convert TT seconds to UT1 days since J2000
    ut1_d = (tt_s + offset_tt2ut1(tt_s))/Tempo.DAY2SEC

    # Retrieve polar motion rotation matrix
    Wt = iers_rot3_itrf_to_tirf(tt_s, m)

    # Compute Earth rotation matrix
    Rt = iers_era_rotm(m, ut1_d)'

    # Assemble rotation
    return Rt*Wt

end



# Miscellaneous Rotations 
# ==============================================

# Low-level function to avoid repeating the same data
function mod_to_gtod3(ttc::Number, m::IERSModel, δΔψ::Number, δΔϵ::Number)

    # Rather than calling the `iers_nutation` function, we manually assembly 
    # the nutation matrix so that we don't have to compute 
    # the nutation components twice (for N and for GAST)

    # Compute mean obliquity at epoch 
    ϵₐ = iers_obliquity(m, ttc)

    # Compute nutation in longitude and obliquity 
    Δψ, Δϵ = iers_nutation_comp(m, ttc)

    # Compute nutation matrix with EOP corrections
    N = angle_to_dcm(ϵₐ, - (Δψ + δΔψ), - (ϵₐ + Δϵ + δΔϵ), :XZX)

    # Compute Greenwich Mean Sidereal Time 
    gmst = iers_gmst(m, ttc)

    # Manually compute the equation of the equinoxes
    gast = gmst + (Δψ + δΔψ)*cos(ϵₐ) + eeq_complementary(m, ttc)

    # Compute sidereal rotation matrix
    R = angle_to_dcm(gast, :Z)

    # Assemble final rotation matrix 
    return R*N

end