
export iers_cip_motion, cip_xy, cip_xys, cip_vector

"""
    iers_cip_motion(m::IERSModel, t::Number, δX::Number=0, δY::Number=0)

Compute the GCRF-to-CIRF rotation matrix, following the IERS Conventions `m`, at time `t` 
expressed in `TT` Julian centuries since J2000. Optional IERS EOP corrections for free-core 
nutation and time dependent effects can be provided via `δX` and `δY`

### References 
- IERS Technical Note No. [36](https://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn36.html) 

### See also 
See also [`cip_xy`](@ref) and [`cip_xys`](@ref).
"""
function iers_cip_motion(m::IERSModel, t::Number, δX::Number=0, δY::Number=0)

    # Compute CIP coordinates 
    X, Y, s = cip_xys(m, t, δX, δY)

    # Form intermediate-to-celestial matrix 
    return xys2m(X, Y, s)

end


"""
    cip_xy(m::IERSModel, t::Number)

Compute the CIP X and Y coordinates, in radians, following the IERS Conventions `m`, at 
time `t`, expressed in `TT` Julian centuries since J2000.

!!! warning 
    The computation of the free-core nutation and time dependent effects are excluded from 
    this model. To achieve the < 1μas accuracy with the IAU 2006/2000 A precession-nutation 
    models, such effects must be included a-posteriori (through δX and δY) using the IERS 
    EOP data.

### References 
- Capitaine N. and Wallace P. T. (2008), Concise CIO based precession-nutation formulations.
- IERS Technical Note No. [21](https://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn21.html)
- IERS Technical Note No. [32](https://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn32.html)
- IERS Technical Note No. [36](https://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn36.html) 

### See also 
See also [`iers_cip_motion`](@ref) and [`cip_xys`](@ref).
"""
cip_xy


"""
    cip_xys(m::IERSModel, t::Number, δX::Number=0, δY::Number=0)

Compute the CIP X, Y and CIO locator s coordinates, in radians, following the IERS 
conventions `m` at time `t`, expressed in `TT` Julian centuries since J2000. Optional EOP 
nutation corrections can be provided via the `δX` and `δY` parameters.

!!! note 
    Because of the small values of the CIP corrections, the CIO locator holds pretty much 
    irrespective on them. Indeed, some reports compute `s` with the corrected CIP 
    coordinates, some do not.

### References 
- Capitaine N. and Wallace P. T. (2008), Concise CIO based precession-nutation formulations
- IERS Technical Note No. [21](https://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn21.html)
- IERS Technical Note No. [32](https://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn32.html)
- IERS Technical Note No. [36](https://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn36.html) 

### See also 
See also [`iers_cip_motion`](@ref) and [`cip_xy`](@ref).
"""
function cip_xys(m::IERSModel, t::Number, δX::Number=0, δY::Number=0)

    # Compute CIP coordinates 
    x, y = cip_xy(m, t)

    # Add EOP corrections
    X = x + δX 
    Y = y + δY

    # Compute CIO locator 
    s = cio_locator(m, t, X, Y)

    return X, Y, s

end


"""
    cip_vector(m::IERSModel, t::Number)

Compute the Celestial Intermediate Pole (CIP) vector, following the IERS Conventions `m` at 
time `t`, expressed in `TT` Julian centuries since J2000.

### References 
- Wallace P. T. and Capitaine N. (2006), Precession-nutation procedures consistent with 
  IAU 2006 resolutions, [DOI: 10.1051/0004-6361:20065897](https://www.aanda.org/articles/aa/abs/2006/45/aa5897-06/aa5897-06.html)

### See also 
See also [`cip_xy`](@ref).
"""
function cip_vector(m::IERSModel, t::Number)
    xs, ys = cip_xy(m, t)
    return SA[xs, ys, sqrt(1 - xs^2 - ys^2)]
end


"""
    cio_locator(m::IERSModel, t::Number, x::Number, y::Number)

Compute the CIO locator `s`, in radians, following the IERS Conventions `m` at time `t`, 
expressed in `TT` Julian centuries since J2000, given the CIP coordinates `x` and `y`.

### References 
- Capitaine N. and Wallace P. T. (2008), Concise CIO based precession-nutation formulations
- IERS Technical Note No. [21](https://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn21.html)
- IERS Technical Note No. [32](https://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn32.html)
- IERS Technical Note No. [36](https://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn36.html) 

### See also 
See also [`iers_cip_motion`](@ref), [`cip_xy`](@ref) and [`cip_xys`](@ref).
"""
function cio_locator(m::IERSModel, t, x, y)
    _cio_locator(m, t, DelaunayArgs(m, t), PlanetaryArgs(m, t)) - x*y/2
end


# 1996 CONVENTIONS
# ============================
function cip_xy(m::IERS1996, t::Number)

    # Compute Delaunay's arguments
    d = DelaunayArgs(m, t)

    # Evaluate X, Y harmonic series 
    xh, yh = _cip_xy(m, t, d)

    # Compute polynomial parts 
    xp = arcsec2rad(@evalpoly(t, 0, 2004.3109, -0.42665, -0.198656, 0.000014))
    yp = arcsec2rad(@evalpoly(t, -0.00013, 0, -22.40992, 0.001836, 0.0011130))

    # Compute mean obliquity at reference epoch 
    ϵ₀ = iers_obliquity(m, 0)

    sΩ, cΩ = sincos(d.Ω) 
    sA, cA = sincos(2*(d.F - d.D + d.Ω))

    # Compute final CEP coordinates 
    x = xp + sin(ϵ₀)*xh + arcsec2rad(0.00204*sΩ + 0.00016*sA)*t^2
    y = yp + yh + arcsec2rad(-0.00231*cΩ - 0.00014*cA)*t^2

    return x, y 

end

function cio_locator(m::IERS1996, t::Number, x::Number, y::Number)

    d = DelaunayArgs(m, t)

    sΩ = sin(d.Ω) 
    sA = sin(2*(d.F - d.D + d.Ω))

    # Compute the polynomial and harmonic parts 
    sp = @evalpoly(t, 0, 0.00385, -0.07259)
    sh = -0.00264*sΩ - 0.00006*sin(2*d.Ω) + t^2*(0.00074*sΩ + 0.00006*sA)

    return arcsec2rad(sp + sh) - x*y/2

end

include("constants/cip1980.jl")
build_series(
    :_cip_xy, :IERS1996, [COEFFS_CIP1980_X, COEFFS_CIP1980_Y]; enable_pargs=false,
    unit_factor=arcsec2rad(1e-4)
)


# 2003 CONVENTIONS
# ============================

function cip_xy(m::IERS2003, t::Number)
    # Extracted from the IAU-2000 bias-precession-nutation matrix 
    return npb2xy(iers_npb(m, t))
end

include("constants/cio2000.jl")
build_series(:_cio_locator, :IERS2003, [COEFFS_CIO2000_S], [COEFFS_CIO2000_SP])


# 2010 CONVENTIONS
# ============================

function cip_xy(m::IERS2010, t::Number)
    
    # Compute Fukushima-Williams angles 
    γ, ϕ, ψ, ϵ = fw_angles(m, t)

    # Compute the IAU 2000 nutation components 
    Δψ, Δϵ = iers_nutation_comp(m, t)

    # Retrieve the CIP coordinates by applying IAU-2006 compatible nutations 
    return fw2xy(γ, ϕ, ψ + Δψ, ϵ + Δϵ)

end

function cip_xy(m::CPNC, t::Number)
    return _cip_xy(m, t, DelaunayArgs(m, t))
end

function cip_xy(::CPND, t::Number)

    μas2rad = 1e-6 * π / 648000

    # Approximated fundamental arguments as linear function of time 
    Ω = 2.182439196616 - 33.7570459536t
    A = -2.776244621014 + 1256.6639307381t

    sΩ, cΩ = sincos(Ω)
    sA, cA = sincos(A)

    x = @evalpoly(t, 0, 2004191898, -429782.9, -198618.34)
    x -= 6844318*sΩ + 523908*sA

    y = @evalpoly(t, 0, 0, -22407275)
    y += 9205236*cΩ + 573033*cA

    x *= μas2rad
    y *= μas2rad

    return x, y
end

function cio_locator(::CPNC, t::Number, x::Number, y::Number)

    # Simplified model! 
    Ω = 2.182439196616 - 33.7570459536t

    s = @evalpoly(t, 0, 3809, 0, -72574.0)
    s -= 2641 * sin(Ω)

    # Transform s from μas to radians 
    s *= 1e-6 * π / 648000

    return s - x * y / 2

end

cio_locator(m::CPND, ::Number, ::Number, ::Number) = 0

include("constants/cip2006c.jl")
build_series(
    :_cip_xy, :CPNC, [COEFFS_CPNC_X, COEFFS_CPNC_Y], [COEFFS_CPNC_XP, COEFFS_CPNC_YP]; 
    enable_pargs=false
)

include("constants/cio2006.jl")
build_series(:_cio_locator, :IERS2010, [COEFFS_CIO2006_S], [COEFFS_CIO2006_SP])


# MISCELLANEOUS
# ============================

""" 
    xys2m(x::Number, y::Number, s::Number)

Compute the GCRF-to-CIRF matrix given the CIP `x`, `y' coordinates and the CIO 
locator `s`, all in radians.

### References
- Wallace P. T. and Capitaine N. (2006), Precession-nutation procedures consistent with 
  IAU 2006 resolutions, [DOI: 10.1051/0004-6361:20065897](https://www.aanda.org/articles/aa/abs/2006/45/aa5897-06/aa5897-06.html) 
- [ERFA](https://github.com/liberfa/erfa/blob/master/src/c2ixys.c) library

"""
function xys2m(x::Number, y::Number, s::Number)
    
    # Retrieves spherical angles E and d. 
    r2 = x^2 + y^2
    E = (r2 > 0) ? atan(y, x) : 0
    d = atan(sqrt(r2 / (1 - r2)))

    # This formulation (compared to simplified versions) ensures 
    # that the resulting matrix is orthonormal.
    return angle_to_dcm(E, d, - (E + s), :ZYZ)

end


"""
    npb2xy(A::AbstractMatrix)

Retrieve the CIP X and Y coordinates from the nutation-precession-bias matrix, in radians.

### References 
- Wallace P. T. and Capitaine N. (2006), Precession-nutation procedures consistent with 
  IAU 2006 resolutions, [DOI: 10.1051/0004-6361:20065897](https://www.aanda.org/articles/aa/abs/2006/45/aa5897-06/aa5897-06.html) 
"""
npb2xy(A::AbstractMatrix) = @inbounds A[3, 1], A[3, 2]