
export orient_era, orient_era_rotm

""" 
    orient_era(m::IERSConventions, t::Number)

Compute the Earth Rotation Angle (ERA), in radians, at time `t` expressed as UT1 days 
since `J2000`, according to the IERS convention `m`.

!!! note 
    In the IERS 1996 conventions, θ is referred to as the Stellar Angle.

### References 
- IERS Technical Note No. [21](https://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn21.html)
- IERS Technical Note No. [36](https://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn36.html) 

### See also 
See also [`orient_era_rotm`](@ref).
"""
function orient_era(::IERSConventions, t::Number)

    # The function uses the fractional UT1 date to gain additional orient_bias_precession_nutation
    # in the computations
    return mod2pi(2π * (mod(t, 1) + 0.7790572732640 + 0.00273781191135448t))

end


"""
    orient_era_rotm(m::IERSConventions, t::Number)

Compute the CIRF-to-TIRF rotation matrix, according to the IERS conventions `m`, at time 
`t` expressed in UT1 days since J2000 

### References
- IERS Technical Note No. [36](https://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn36.html) 

### See also 
See also [`orient_era`](@ref).
"""
function orient_era_rotm(m::IERSConventions, t::Number) 
    return angle_to_dcm(orient_era(m, t), :Z)
end