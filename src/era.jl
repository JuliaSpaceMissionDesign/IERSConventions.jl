
export iers_era, iers_era_rotm

""" 
    iers_era(m::IERSModel, ut1_d::Number)

Compute the Earth Rotation Angle (ERA), in radians, at time `ut1_d` expressed as UT1 days 
since `J2000`, according to the IERS convention `m`.

!!! note 
    In the IERS 1996 conventions, θ is referred to as the Stellar Angle.

### References 
- IERS Technical Note No. [21](https://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn21.html)
- IERS Technical Note No. [36](https://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn36.html) 

### See also 
See also [`iers_era_rotm`](@ref).
"""
function iers_era(::IERSModel, ut1_d::Number)
    # The function uses the fractional UT1 date to gain additional precision
    # in the computations

    # Uses integer divsion because the `mod` and `rem` functions don't work with ForwardDiff
    f = ut1_d - (ut1_d ÷ 1)
    return rem2pi(2π * (f + 0.7790572732640 + 0.00273781191135448ut1_d), RoundDown)

end


"""
    iers_era_rotm(m::IERSModel, ut1_d::Number)

Compute the CIRF-to-TIRF rotation matrix, according to the IERS conventions `m`, at time 
`ut1_d` expressed in UT1 days since J2000 

### References
- IERS Technical Note No. [36](https://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn36.html) 

### See also 
See also [`iers_era`](@ref).
"""
function iers_era_rotm(m::IERSModel, ut1_d::Number) 
    return angle_to_dcm(iers_era(m, ut1_d), :Z)
end
