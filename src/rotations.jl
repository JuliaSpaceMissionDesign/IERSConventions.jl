
# Equinox-based transformations from GCRF
# ==============================================

"""
    iers_rot3_gcrf_to_mod(t::Number, m::IERSModel=iers2010b)

Compute the rotation matrix from the Geocentric Celestial Reference Frame (GCRF) to 
the Mean-of-Date (MOD) at time `t`, expressed in TT seconds since `J2000`.

!!! note 
    The Mean-of-Date axes are obtained by applying the frame bias and precession matrix. 
    For this reason, if the [`iers1996`](@ref) conventions are used, the rotation is 
    actually computed starting from the MEME2000 rather than the GCRF.  

### See also 
See also [`iers_pb`](@ref) and [`iers_rot3_itrf_to_mod`](@ref).

"""
function iers_rot3_gcrf_to_mod(t::Number, m::IERSModel=iers2010b)
    return iers_pb(m, t/Tempo.CENTURY2SEC)
end


"""
    iers_rot3_gcrf_to_tod(t::Number, m::IERSModel=iers2010b)

Compute the rotation matrix from the Geocentric Celestial Reference Frame (GCRF) to 
the True-of-Date (TOD) at time `t`, expressed in TT seconds since `J2000`.

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
function iers_rot3_gcrf_to_tod(t::Number, m::IERSModel=iers2010b)

    # Retrieve the EOP corrections to the nutation in longitude and obliquity
    δΔψ = eop_δΔψ(m, t)
    δΔϵ = eop_δΔϵ(m, t)

    # Compute Nutation-Precession-Bias matrix
    return iers_npb(m, t/Tempo.CENTURY2SEC, δΔψ, δΔϵ)

end


"""
    iers_rot3_gcrf_to_gtod(t::Number, m::IERSModel=iers2010b)

Compute the rotation matrix from the Geocentric Celestial Reference Frame (GCRF) to 
the Greenwich True-of-Date (GTOD) at time `t`, expressed in TT seconds since `J2000`.

!!! note 
    If the [`iers1996`](@ref) conventions are used, the rotation is actually computed 
    starting from the MEME2000 rather than the GCRF.  

!!! note 
    The EOP nutation corrections are only used for the [`iers2003a`](@ref) and 
    [`iers2010a`](@ref) models. 

## See also 
See also [`mod_to_gtod3`](@ref) and [`iers_rot3_itrf_to_gtod`](@ref).
"""
function iers_rot3_gcrf_to_gtod(t::Number, m::IERSModel=iers2010b)

    # Retrieve the EOP corrections to the nutation in longitude and obliquity
    δΔψ = eop_δΔψ(m, t)
    δΔϵ = eop_δΔϵ(m, t)

    # Retrieve precession-bias matrix
    PB = iers_pb(m, tt_c)

    # Retrieve MOD-to-GTOD rotation 
    RN = mod_to_gtod3(t, m, δΔψ, δΔϵ)
    return RN*PB

end


"""
    iers_rot3_gcrf_to_pef(t::Number, m::IERSModel=iers2010b)

Compute the rotation matrix from the Geocentric Celestial Reference Frame (GCRF) to 
the Pseudo-Earth Fixed (PEF) at time `t`, expressed in TT seconds since `J2000`.

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
function iers_rot3_gcrf_to_pef(t::Number, m::IERSModel=iers2010b)

    # Compute the GCRF to GTOD rotation matrix
    RNPB = iers_rot3_gcrf_to_gtod(t, m)

    # Compute the TIO locator matrix (from GTOD to PEF)
    sp = tio_locator(m, t/Tempo.CENTURY2SEC)

    # Assemble the global rotation
    return angle_to_dcm(sp, :Z)*RNPB

end

# TODO: add TEME rotations


# Equinox-based transformations from ITRF
# ==============================================

"""
    iers_rot3_itrf_to_pef(t::Number, m::IERSModel=iers2010b)

Compute the rotation matrix from the International Terrestrial Reference Frame (ITRF) to 
the Pseudo-Earth Fixed (PEF) at time `t`, expressed in TT seconds since `J2000`.

!!! note 
    The [`CPNd`](@ref) model returns an identity rotation because it neglects the 
    effects of polar motion. # TODO: specify magnitude!

## See also 
See also [`iers_rot3_gcrf_to_pef`](@ref).
"""
function iers_rot3_itrf_to_pef(t::Number, ::IERSModel=iers2010b)

    # Retrieve pole coordinates 
    xₚ = eop_xp(m, t)
    yₚ = eop_yp(m, t)

    # Compute the polar motion rotation without the TIO contribution 
    # Note: this equation is wrongly reported on Vallado Fundamentals of Astrodyn. 4th Ed. 
    return angle_to_dcm(yₚ, xₚ, :XY)

end

iers_rot3_itrf_to_pef(::Number, ::CPND) = DCM(1I)


"""
    iers_rot3_itrf_to_gtod(t::Number, m::IERSModel=iers2010b)

Compute the rotation matrix from the International Terrestrial Reference Frame (ITRF) to 
the Greenwich True-of-Date (GTOD) at time `t`, expressed in TT seconds since `J2000`.

!!! note 
    The [`CPNd`](@ref) model returns an identity rotation because it neglects the 
    effects of polar motion. # TODO: specify magnitude!

!!! note 
    For the [`iers1996`](@ref) model, the GTOD and PEF axes are equal because the IERS 
    1996 conventions do not account for the TIO locator effects.

## See also 
See also [`iers_rot3_itrf_to_pef`](@ref) and [`iers_rot3_gcrf_to_gtod`](@ref).
"""
function iers_rot3_itrf_to_gtod(t::Number, m::IERSModel=iers2010b)

    # Compute ITRF to PEF rotation
    W  = iers_rot3_itrf_to_pef(t, m)

    # Compute TIO locator rotation 
    sp = tio_locator(m, t/Tempo.CENTURY2SEC)
    return angle_to_dcm(-sp, :Z)*W

end


"""
    iers_rot3_itrf_to_tod(t::Number, m::IERSModel=iers2010b)

Compute the rotation matrix from the International Terrestrial Reference Frame (ITRF) to 
the True-of-Date (TOD) at time `t`, expressed in TT seconds since `J2000`.

!!! note 
    The EOP corrections in longitude are only used for the [`iers2003a`](@ref) and 
    [`iers2010a`](@ref) models. 

## See also 
See also [`iers_rot3_itrf_to_gtod`](@ref) and [`iers_rot3_gcrf_to_tod`](@ref).
"""
function iers_rot3_itrf_to_tod(t::Number, m::IERSModel=iers2010b)

    # Compute EOP nutation corrections 
    δΔψ = eop_δΔψ(m, t)

    # Compute ITRF to GTOD rotation 
    W = iers_rot3_itrf_to_gtod(m, t)

    # Compute GTOD to TOD rotation matrix 
    R = angle_to_dcm(-iers_gast(m, t/Tempo.CENTURY2SEC, δΔψ), :Z)
    return R*W

end


"""
    iers_rot3_itrf_to_mod(t::Number, m::IERSModel=iers2010b)

Compute the rotation matrix from the International Terrestrial Reference Frame (ITRF) to 
the Mean-of-Date (MOD) at time `t`, expressed in TT seconds since `J2000`.

!!! note 
    The EOP nutation corrections are only used for the [`iers2003a`](@ref) and 
    [`iers2010a`](@ref) models. 

## See also 
See also [`iers_rot3_itrf_to_tod`](@ref) and [`iers_rot3_gcrf_to_mod`](@ref).
"""
function iers_rot3_itrf_to_mod(t::Number, m::IERSModel=iers2010b)

    # Compute ITRF to GTOD rotation 
    W = iers_rot3_itrf_to_gtod(m, t)

    # Retrieve the EOP corrections to the nutation in longitude and obliquity
    δΔψ = eop_δΔψ(m, t)
    δΔϵ = eop_δΔϵ(m, t)

    # Compute GTOD-to-MOD rotation
    RNt = mod_to_gtod3(t, m, δΔψ, δΔϵ)' 
    return RNt*W 

end


# CIO-based transformations from GCRF 
# ==============================================

"""
    iers_rot3_gcrf_to_cirf(t::Number, m::IERSModel=iers2010b)

Compute the rotation matrix from the Geocentric Celestial Reference Frame (GCRF) to 
the Celestial Intermediate Reference Frame (CIRF) at time `t`, expressed in TT seconds 
since `J2000`, following the IERS Conventions `m`. 

!!! note 
    EOP corrections to the CIP coordinates (δX, δY) are only added in the [`iers2003a`](@ref)
    and [`iers2010a`](@ref) models.

### References 
- IERS Technical Note No. [36](https://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn36.html) 

### See also 
See also [`iers_rot3_gcrf_to_tirf`](@ref) and [`iers_rot3_gcrf_to_itrf`](@ref).
"""
function iers_rot3_gcrf_to_cirf(t::Number, m::IERSModel=iers2010b)
    return iers_cip_motion(m, t/Tempo.CENTURY2SEC, eop_δX(m, t), eop_δY(m, t))
end

"""
    iers_rot3_gcrf_to_tirf(t::Number, m::IERSModel=iers2010b)

Compute the rotation matrix from the Geocentric Celestial Reference Frame (GCRF) to 
the Terrestrial Intermediate Reference Frame (TIRF) at time `t`, expressed in TT seconds 
since `J2000`, following the IERS Conventions `m`. 

!!! note 
    EOP corrections to the CIP coordinates (δX, δY) are only added in the [`iers2003a`](@ref)
    and [`iers2010a`](@ref) models.

### References 
- IERS Technical Note No. [36](https://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn36.html) 

### See also 
See also [`iers_rot3_gcrf_to_cirf`](@ref) and [`iers_rot3_gcrf_to_itrf`](@ref).
"""
function iers_rot3_gcrf_to_tirf(t::Number, m::IERSModel=iers2010b)

    # Convert TT seconds to UT1 days since J2000
    ut1_d = (t + offset_tt2ut1(t))/Tempo.DAY2SEC

    # Compute GCRF to CIRF rotation (CIP motion)
    Q = iers_rot3_gcrf_to_cirf(t, m)

    # Compute the ERA rotation matrixes
    R = iers_era_rotm(m, ut1_d)

    return R*Q

end


"""
    iers_rot3_gcrf_to_itrf(t::Number, m::IERSModel=iers2010b)

Compute the rotation matrix from the Geocentric Celestial Reference Frame (GCRF) to 
the International Terrestrial Reference Frame (ITRF) at time `t`, expressed in TT seconds 
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
function iers_rot3_gcrf_to_itrf(t::Number, m::IERSModel=iers2010b)

    # Compute the GCRF to TIRF rotation matrix
    RQ = iers_rot3_gcrf_to_tirf(t, m)

    # Retrieve pole coordinates 
    xₚ = eop_xp(m, t)
    yₚ = eop_yp(m, t)

    # Compute the polar motion rotation
    W = iers_polar_motion(m, xₚ, yₚ, t/Tempo.CENTURY2SEC)

    return W*RQ

end


# CIO-based transformations from ITRF 
# ==============================================

"""
    iers_rot3_itrf_to_tirf(t::Number, m::IERSModel=iers2010b)

Compute the rotation matrix from the International Terrestrial Reference Frame (ITRF) to 
the Terrestrial Intermediate Reference Frame (TIRF) at time `t`, expressed in TT seconds 
since `J2000`, following the IERS Conventions `m`. 

!!! note 
    Polar motion is neglected in the [`CPNd`](@ref) model. 

### References 
- IERS Technical Note No. [36](https://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn36.html) 

### See also 
See also [`iers_rot3_gcrf_to_tirf`](@ref) and [`iers_rot3_itrf_to_cirf`](@ref).
"""
function iers_rot3_itrf_to_tirf(t::Number, m::IERSModel=iers2010b)

    # Retrieve pole coordinates 
    xₚ = eop_xp(m, t)
    yₚ = eop_yp(m, t)

    # Compute the polar motion rotation
    return iers_polar_motion(m, xₚ, yₚ, t/Tempo.CENTURY2SEC)'

end

iers_rot3_itrf_to_tirf(::Number, ::CPND) = DCM(1I)


"""
    iers_rot3_itrf_to_cirf(t::Number, m::IERSModel=iers2010b)

Compute the rotation matrix from the International Terrestrial Reference Frame (ITRF) to 
the Celestial Intermediate Reference Frame (CIRF) at time `t`, expressed in TT seconds 
since `J2000`, following the IERS Conventions `m`. 

!!! note 
    Polar motion is neglected in the [`CPNd`](@ref) model. 

### References 
- IERS Technical Note No. [36](https://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn36.html) 

### See also 
See also [`iers_rot3_gcrf_to_cirf`](@ref) and [`iers_rot3_itrf_to_tirf`](@ref).
"""
function iers_rot3_itrf_to_cirf(t::Number, m::IERSModel=iers2010b)

    # Convert TT seconds to UT1 days since J2000
    ut1_d = (t + offset_tt2ut1(t))/Tempo.DAY2SEC

    # Retrieve polar motion rotation matrix
    Wt = iers_rot3_itrf_to_tirf(t, m)

    # Compute Earth rotation matrix
    Rt = iers_era_rotm(m, ut1_d)'

    # Assemble rotation
    return Rt*Wt

end



# Miscellaneous Rotations 
# ==============================================

# Low-level function to avoid repeating the same data
function mod_to_gtod3(t::Number, m::IERSModel, δΔψ::Number, δΔϵ::Number)

    # Rather than calling the `iers_nutation` function, we manually assembly 
    # the nutation matrix so that we don't have to compute 
    # the nutation components twice (for N and for GAST)

    # Convert TT seconds to TT Julian centuries since J2000
    tt_c = t/Tempo.CENTURY2SEC

    # Compute mean obliquity at epoch 
    ϵₐ = iers_obliquity(m, tt_c)

    # Compute nutation in longitude and obliquity 
    Δψ, Δϵ = iers_nutation_comp(m, tt_c)

    # Compute nutation matrix with EOP corrections
    N = angle_to_dcm(ϵₐ, - (Δψ + δΔψ), - (ϵₐ + Δϵ + δΔϵ), :XZX)

    # Compute Greenwich Mean Sidereal Time 
    gmst = iers_gmst(m, tt_c)

    # Manually compute the equation of the equinoxes
    gast = gmst + (Δψ + δΔψ)*cos(ϵₐ) + eeq_complementary(m, tt_c)

    # Compute sidereal rotation matrix
    R = angle_to_dcm(gast, :Z)

    # Assemble final rotation matrix 
    return R*N

end