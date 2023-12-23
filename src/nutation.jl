
export iers_nutation, iers_nutation_comp


"""
    iers_nutation(m::IERSConventions, t::Number, δΔψ::Number=0, δΔϵ::Number=0)

Compute the nutation matrix that rotates a vector from Mean-of-Date (MOD) to True-of-Date 
(TOD) axes following the IERS convention `m`, at time `t` expressed in `TT` Julian 
Centuries since `J2000`. 

Optional EOP nutation corrections can be provided via the `δΔψ` and `δΔϵ` parameters.

### References 
- Wallace P. T. and Capitaine N. (2006), Precession-nutation procedures consistent with 
  IAU 2006 resolutions, [DOI: 10.1051/0004-6361:20065897](https://www.aanda.org/articles/aa/abs/2006/45/aa5897-06/aa5897-06.html) 

### See also 
See also [`iers_nutation_comp`](@ref) and [`iers_obliquity`](@ref). 
"""
function iers_nutation(m::IERSConventions, t::Number, δΔψ::Number=0, δΔϵ::Number=0)
    
    # Compute mean obliquity at epoch 
    ϵₐ = iers_obliquity(m, t)

    # Compute nutation in longitude and obliquity 
    Δψ, Δϵ = iers_nutation_comp(m, t)

    # Compute nutation matrix with EOP corrections
    return angle_to_dcm(ϵₐ, - (Δψ + δΔψ), - (ϵₐ + Δϵ + δΔϵ), :XZX)
    
end


"""
    iers_nutation_comp(m::IERSConventions, t::Number)

Compute the nutation components in longitude and obliquity for the IERS convention `m`, in 
radians, at time `t` expressed in `TT` Julian Centuries since `J2000`.

!!! note 
    For the **IAU 2006A** model, the function strictly follows the SOFA implementation. It 
    first computes the IAU 2000A nutation, then applies adjustments for the consequences of 
    the change in obliquity from the IAU 1980 ecliptic to the IAU 2006 ecliptic and (ii) 
    for the secular variation in the Earth's dynamical form factor J2. These corrections 
    ensure that the IAU 2000A nutation is consistent with the IAU 2006 precession model. 
    Please note that the coefficients available on the IERS tables already include those 
    corrections, and are retrieved by multiplying the amplitudes of the SOFA nutation in 
    longitude coefficients by 1.00000047. 

!!! note 
    The expressions of these components for the [`CPNc`](@ref) and [`CPNd`](@ref) models are 
    indirectly computed from their CIP series expansion.

!!! warning 
    The computation of the free-core nutation and time dependent effects are excluded from 
    this model. To achieve the < 1μas accuracy with the IAU 2006/2000 A precession-nutation 
    models, such effects must be included a-posteriori (through δΔψ and δΔϵ) using the IERS 
    EOP data.

### References 
- IERS Technical Note No. [21](https://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn21.html)
- IERS Technical Note No. [32](https://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn32.html)
- IERS Technical Note No. [36](https://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn36.html) 
- Wallace P. T. and Capitaine N. (2006), Precession-nutation procedures consistent with 
  IAU 2006 resolutions, [DOI: 10.1051/0004-6361:20065897](https://www.aanda.org/articles/aa/abs/2006/45/aa5897-06/aa5897-06.html)
- Capitaine N. and Wallace P. T. (2008), Concise CIO based precession-nutation formulations
- ERFA [nut06a](https://github.com/liberfa/erfa/blob/master/src/nut06a.c) and 
  [nut00b](https://github.com/liberfa/erfa/blob/master/src/nut00b.c) functions 

### See also 
See also [`iers_nutation`](@ref)
"""
iers_nutation_comp


# 1996 CONVENTIONS
# ============================

function iers_nutation_comp(m::IERS1996, t::Number)
    
    # Compute Delaunay's arguments 
    dargs = DelaunayArgs(m, t)

    # Compute the nutation components in longitude and obliquity
    Δψ, Δϵ = _nut_components(m, t, dargs)

end

include("constants/nut1980.jl")

build_series(
    :_nut_components, :IERS1996, [COEFFS_NUT1980_ψ, COEFFS_NUT1980_ϵ]; 
    enable_pargs=false, unit_factor=arcsec2rad(1e-4)
);


# 2003 CONVENTIONS
# ============================

function iers_nutation_comp(m::IERS2003A, t::Number)

    # Compute Delaunay's arguments 
    dargs = DelaunayArgs(m, t)

    # Compute Planetary's arguments 
    pargs = PlanetaryArgs(m, t)

    # Compute the nutation components in longitude and obliquity
    Δψ, Δϵ = _nut_components(m, t, dargs, pargs)
    return Δψ, Δϵ

end

function iers_nutation_comp(m::IERS2003B, t::Number)

    # Compute Delaunay's arguments 
    dargs = DelaunayArgs(m, t)

    # Compute the nutation components
    Δψ_ls, Δϵ_ls = _nut_components(m, t, dargs) 

    # Add the offeset to account for the truncated planetary contributions 
    δψ_pl = arcsec2rad(-0.135 * 1e-3)
    δϵ_pl = arcsec2rad(0.388 * 1e-3)

    # Return the nutation components in longitude and obliquity
    return Δψ_ls + δψ_pl, Δϵ_ls + δϵ_pl

end

include("constants/nut2000a.jl")
include("constants/nut2000b.jl")

build_series(:_nut_components, :IERS2003A, [COEFFS_NUT2000A_ψ, COEFFS_NUT2000A_ϵ]);
build_series(
    :_nut_components, :IERS2003B, [COEFFS_NUT2000B_ψ, COEFFS_NUT2000B_ϵ]; 
    enable_pargs=false
);


# 2010 CONVENTIONS
# ============================

function iers_nutation_comp(::IERS2010A, t::Number)

    # Computes IAU 2000A nutation components from luni-solar 
    # and planetary terms of the Mathews et al. (2002) series
    Δψₐ, Δϵₐ = iers_nutation_comp(iers2003a, t)

    # Factor correcting the secular variation of J2 
    fj2 = -2.7774e-6 * t # t = TT 

    # Applies P03 Nutation Corrections from WC06 (2006)
    Δψ = (1 + fj2 + 0.4697e-6) * Δψₐ
    Δϵ = (1 + fj2) * Δϵₐ
    
    # Return the nutation components in longitude and obliquity
    return Δψ, Δϵ

end

iers_nutation_comp(::IERS2010B, t) = iers_nutation_comp(iers2003b, t)

function iers_nutation_comp(m::CPNC, t::Number)

    # Compute CIP coordinates
    x, y = cip_xy(m, t)
    
    # Compute mean obliquity at epoch 
    ϵₐ = iers_obliquity(m, t)

    # Compute the precession angles 
    ϵ₀, ψₐ, ωₐ, _ = precession_angles_rot4(m, t)
    se, ce = sincos(ϵ₀)

    # Frame bias quantities
    η₀ = arcsec2rad(-0.0068192)
    ξ₀ = arcsec2rad(-0.016617)

    b = se*ce*ψₐ

    a = y - η₀ - (ωₐ - ϵ₀) + se*ce/2*ψₐ^2
    c = x - ξ₀ - ψₐ*se + se/6*ψₐ^3 - ψₐ*ce*a 

    # Compute nutation in longitude and obliquity
    Δψ = c/(se + (ϵₐ - ϵ₀)*ce + ψₐ*ce*b)
    Δϵ = a + b*Δψ

    return Δψ, Δϵ

end

function iers_nutation_comp(m::CPND, t::Number)

    # Compute CIP coordinates 
    x, y = cip_xy(m, t)

    # Retrieve precession angles 
    ϵ₀, ψₐ, _, _ = precession_angles_rot4(m, t)
    se, ce = sincos(ϵ₀)

    # Compute nutation in longitude and obliquity
    Δψ = x/se - ψₐ
    Δϵ = y + se*ce/2*ψₐ^2

    return Δψ, Δϵ

end