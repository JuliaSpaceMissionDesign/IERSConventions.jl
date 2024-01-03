
export iers_bias, iers_pb, iers_npb

"""
    iers_bias(m::IERSModel, t::Number)

Compute the frame bias matrix, which transform vectors from the GCRF axes to the Mean 
Equinox and Mean Equator of J2000 (MEME2000) axes. 

!!! note 
    Since in the IERS 1996 conventions the bias matrix was still undefined, the returned 
    matrix is the identity. 

### References 
- IERS Technical Note No. [21](https://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn21.html)
- IERS Technical Note No. [32](https://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn32.html)
- IERS Technical Note No. [36](https://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn36.html) 

### See also 
See also [`iers_precession`](@ref), [`iers_pb`](@ref) and [`iers_npb`](@ref).
"""
iers_bias


"""
    iers_pb(m::IERSModel, t::Number)

Compute the precession-bias (PB) matrix, which transforms vectors from the GCRF axes to
Mean-of-Date (MOD) axes, at time `t` expressed in `TT` Julian centuries since `J2000`, according 
to the IERS convention `m`

### References 
- Wallace P. T. and Capitaine N. (2006), Precession-nutation procedures consistent with 
  IAU 2006 resolutions, [DOI: 10.1051/0004-6361:20065897](https://www.aanda.org/articles/aa/abs/2006/45/aa5897-06/aa5897-06.html) 
- IERS Technical Note No. [21](https://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn21.html)
- IERS Technical Note No. [32](https://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn32.html)
- IERS Technical Note No. [36](https://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn36.html) 

### See also 
See also [`iers_bias`](@ref), [`iers_precession`](@ref) and [`iers_npb`](@ref).
"""
iers_pb


"""
    iers_npb(m::IERSModel, t::Number, δΔψ=0, δΔϵ=0)

Compute the nutation-bias-precession (NPB) matrix, which transforms vectors from the GCRF 
to True-of-Date (TOD) axes, at time `t` expressed in `TT` Julian centuries since `J2000`, 
according to the IERS convention `m`

### References 
- Wallace P. T. and Capitaine N. (2006), Precession-nutation procedures consistent with 
  IAU 2006 resolutions, [DOI: 10.1051/0004-6361:20065897](https://www.aanda.org/articles/aa/abs/2006/45/aa5897-06/aa5897-06.html) 
- IERS Technical Note No. [21](https://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn21.html)
- IERS Technical Note No. [32](https://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn32.html)
- IERS Technical Note No. [36](https://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn36.html) 

### See also 
See also [`iers_nutation`](@ref), [`iers_precession`](@ref) and [`iers_pb`](@ref).
"""
iers_npb


# 1996 CONVENTIONS
# ============================

iers_bias(::IERS1996, ::Number) = DCM(1I)

iers_pb(m::IERS1996, t::Number) = iers_precession(m, t)

function iers_npb(m::IERS1996, t::Number, δΔψ::Number=0, δΔϵ::Number=0)

    # Compute precession and nutation matrices
    rp = iers_precession(m, t)
    rn = iers_nutation(m, t, δΔψ, δΔϵ)

    # Form Precession-Nutation matrix 
    return rn * rp
end


# 2003 CONVENTIONS
# ============================

function iers_bias(m::IERS2003, ::Number)

    # Compute shift in the origin 
    δα₀ = arcsec2rad(-0.0146)

    ϵ₀ = iers_obliquity(m, 0)
    δψ = arcsec2rad(-0.041775)

    # Compute offsets from the ICRS pole
    ξ₀ = δψ * sin(ϵ₀)
    η₀ = arcsec2rad(-0.0068192)

    return angle_to_dcm(δα₀, ξ₀, -η₀, :ZYX)
end

function iers_pb(m::IERS2003, t::Number)

    # Compute bias and precession matrices
    rb = iers_bias(m, t)
    rp = iers_precession(m, t)

    # Form Bias-Precession matrix 
    return rp * rb
end

function iers_npb(m::IERS2003, t::Number, δΔψ::Number=0, δΔϵ::Number=0)

    # Compute bias and precession matrices
    rb = iers_bias(m, t)
    rp = iers_precession(m, t)

    # Compute the nutation matrix with added EOP corrections
    rn = iers_nutation(m, t, δΔψ, δΔϵ)

    # Form Bias-Precession matrix 
    return rn * rp * rb
end


# 2010 CONVENTIONS
# ============================

function iers_bias(m::IERS2010, ::Number)

    # Computes Fukushima-Williams angles at J2000
    γ, ϕ, ψ, ϵ = fw_angles(m, 0)

    # Compute the Bias matrix
    return fw_matrix(γ, ϕ, ψ, ϵ)
end

function iers_pb(m::IERS2010, t::Number)

    # Computes Fukushima-Williams angles at epoch
    γ, ϕ, ψ, ϵ = fw_angles(m, t)

    # Compute the Bias-precession matrix
    return fw_matrix(γ, ϕ, ψ, ϵ)
end

function iers_npb(m::IERS2010, t::Number, δΔψ::Number=0, δΔϵ::Number=0)

    # Computes Fukushima-Williams angles
    γ, ϕ, ψ, ϵ = fw_angles(m, t)

    # Computes IAU 2000 nutation components 
    Δψ, Δϵ = iers_nutation_comp(m, t)

    # Compute the Bias-precession-nutation matrix by applying 
    # the IAU-2006 compatible nutations with added EOP nutation corrections
    return fw_matrix(γ, ϕ, ψ + Δψ + δΔψ, ϵ + Δϵ + δΔϵ)
end
