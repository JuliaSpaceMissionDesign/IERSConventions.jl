
# NOTE: we disable flat extrapolation so we don't have to perform the check 
# on the interpolation boundaries twice.  


"""
    eop_δΔψ(m::IERSModel, tt_c::Number)

Interpolate and retrieve the EOP nutation correction in longitude `δΔψ`, in radians, 
at time `tt_c` expressed in `TT` Julian centuries since J2000 for the IERS convention `m`.

!!! note 
    If time is outside the boundaries for which there is available EOP data, a zero
    value is returned.

### See also 
See also [`eop_δΔϵ`](@ref), [`eop_δX`](@ref) and [`eop_δY`](@ref).
"""
function eop_δΔψ(::IERS2010A, tt_c::Number)
    if eop_check_time(tt_c)        
        return π/648000*interpolate(IERS_EOP.nut2010.δΔψ, tt_c, false)
    else
        return 0.0*zero(tt_c)
    end 
end

function eop_δΔψ(::IERS2003A, tt_c::Number)
    if eop_check_time(tt_c)        
        return π/648000*interpolate(IERS_EOP.nut2003.δΔψ, tt_c, false)
    else
        return 0.0*zero(tt_c)
    end 
end

function eop_δΔψ(::IERS1996, tt_c::Number)
    if eop_check_time(tt_c)        
        return π/648000*interpolate(IERS_EOP.nut1996.δΔψ, tt_c, false)
    else
        return 0.0*zero(tt_c)
    end 
end

eop_δΔψ(::IERSModel, ::Number) = 0


"""
    eop_δΔϵ(m::IERSModel, tt_c::Number)

Interpolate and retrieve the EOP nutation correction in obliquity `δΔϵ`, in radians, 
at time `tt_c` expressed in `TT` Julian centuries since J2000 for the IERS convention `m`.

!!! note 
    If time is outside the boundaries for which there is available EOP data, a zero
    value is returned.

### See also 
See also [`eop_δΔψ`](@ref), [`eop_δX`](@ref) and [`eop_δY`](@ref).
"""
function eop_δΔϵ(::IERS2010A, tt_c::Number)
    if eop_check_time(tt_c)        
        return π/648000*interpolate(IERS_EOP.nut2010.δΔϵ, tt_c, false)
    else
        return 0.0*zero(tt_c)
    end 
end

function eop_δΔϵ(::IERS2003A, tt_c::Number)
    if eop_check_time(tt_c)        
        return π/648000*interpolate(IERS_EOP.nut2003.δΔϵ, tt_c, false)
    else
        return 0.0*zero(tt_c)
    end 
end

function eop_δΔϵ(::IERS1996, tt_c::Number)
    if eop_check_time(tt_c)         
        return π/648000*interpolate(IERS_EOP.nut1996.δΔϵ, tt_c, false)
    else
        return 0.0*zero(tt_c)
    end 
end

eop_δΔϵ(::IERSModel, ::Number) = 0


"""
    eop_δX(m::IERSModel, tt_c::Number)

Interpolate and retrieve the CIP `δX` correction, in radians, at time `tt_c` expressed in 
`TT` Julian centuries since J2000 for the IERS convention `m`.

!!! note 
    If time is outside the boundaries for which there is available EOP data, a zero
    value is returned.

### See also 
See also [`eop_δΔψ`](@ref), [`eop_δΔϵ`](@ref) and [`eop_δY`](@ref).
"""
function eop_δX(::IERS2010A, tt_c::Number)
    if eop_check_time(tt_c)        
        return π/648000*interpolate(IERS_EOP.nut2010.δX, tt_c, false)
    else
        return 0.0*zero(tt_c)
    end 
end

function eop_δX(::IERS2003A, tt_c::Number)
    if eop_check_time(tt_c)        
        return π/648000*interpolate(IERS_EOP.nut2003.δX, tt_c, false)
    else
        return 0.0*zero(tt_c)
    end 
end

function eop_δX(::IERS1996, tt_c::Number)
    if eop_check_time(tt_c)        
        return π/648000*interpolate(IERS_EOP.nut1996.δX, tt_c, false)
    else
        return 0.0*zero(tt_c)
    end 
end

eop_δX(::IERSModel, ::Number) = 0


"""
    eop_δY(m::IERSModel, tt_c::Number)

Interpolate and retrieve the CIP `δY` correction, in radians, at time `tt_c` expressed in 
`TT` Julian centuries since J2000 for the IERS convention `m`.

!!! note 
    If time is outside the boundaries for which there is available EOP data, a zero
    value is returned.

### See also 
See also [`eop_δΔψ`](@ref), [`eop_δΔϵ`](@ref) and [`eop_δX`](@ref).
"""
function eop_δY(::IERS2010A, tt_c::Number)
    if eop_check_time(tt_c)        
        return π/648000*interpolate(IERS_EOP.nut2010.δY, tt_c, false)
    else
        return 0.0*zero(tt_c)
    end 
end

function eop_δY(::IERS2003A, tt_c::Number)
    if eop_check_time(tt_c)        
        return π/648000*interpolate(IERS_EOP.nut2003.δY, tt_c, false)
    else
        return 0.0*zero(tt_c)
    end 
end

function eop_δY(::IERS1996, tt_c::Number)
    if eop_check_time(tt_c)        
        return π/648000*interpolate(IERS_EOP.nut1996.δY, tt_c, false)
    else
        return 0.0*zero(tt_c)
    end 
end

eop_δY(::IERSModel, ::Number) = 0


"""
    eop_xp(m::IERSModel, tt_c::Number)

Interpolate and retrieve the pole `xₚ` coordinate, in radians, at time `tt_c` expressed in 
`TT` Julian centuries since J2000 for the IERS convention `m`.

!!! note 
    If time is outside the boundaries for which there is available EOP data, a zero
    value is returned.

### See also 
See also [`eop_yp`](@ref).
"""
function eop_xp(::IERSModel, tt_c::Number)
    if eop_check_time(tt_c)        
        return π/648000*interpolate(IERS_EOP.xp, tt_c, false)
    else
        return 0.0*zero(tt_c)
    end 
end

eop_xp(::CPND, ::Number) = 0


"""
    eop_yp(m::IERSModel, tt_c::Number)

Interpolate and retrieve the pole `yₚ` coordinate, in radians, at time `tt_c` expressed in 
`TT` Julian centuries since J2000 for the IERS convention `m`.

!!! note 
    If time is outside the boundaries for which there is available EOP data, a zero
    value is returned.

### See also 
See also [`eop_xp`](@ref).
"""
function eop_yp(::IERSModel, tt_c::Number)
    if eop_check_time(tt_c)        
        return π/648000*interpolate(IERS_EOP.yp, tt_c, false)
    else
        return 0.0*zero(tt_c)
    end 
end

eop_yp(::CPND, ::Number) = 0


"""
    offset_tt2ut1(tt_s::Number)

Return the TT-to-UT1 offset, in seconds, at `tt_s` expressed in TT seconds since J2000.

!!! note 
    If time is outside the boundaries for which there is available EOP data, a zero
    value is returned.
"""
function offset_tt2ut1(tt_s)

    tt_c = tt_s/Tempo.CENTURY2SEC

    if eop_check_time(tt_c)        
        return interpolate(IERS_EOP.ut1_tt, tt_c, false)
    else
        return 0.0*zero(tt_c)
    end 

end

# Check whether the input time is within the boundaries of the EOP data
function eop_check_time(tt_c::Number)
    eop_check_init()
    return IERS_EOP_DATA.cent_TT[1] <= tt_c <= IERS_EOP_DATA.cent_TT[end]
end

