
"""
    eop_δΔψ(m::IERSModel, t::Number)

Interpolate and retrieve the EOP nutation correction in longitude `δΔψ`, in radians, 
at time `t` expressed in `TT` Julian centuries since J2000 for the IERS convention `m`.
"""
function eop_δΔψ(::IERS2010A, t::Number)
    eop_check_init()
    return interpolate(IERS_EOP.nut2010.δΔψ, t)
end

function eop_δΔψ(::IERS2003A, t::Number)
    eop_check_init()
    return interpolate(IERS_EOP.nut2003.δΔψ, t)
end

function eop_δΔψ(::IERS1996, t::Number)
    eop_check_init()
    return interpolate(IERS_EOP.nut1996.δΔψ, t)
end

eop_δΔψ(::IERSModel, ::Number) = 0


"""
    eop_δΔϵ(m::IERSModel, t::Number)

Interpolate and retrieve the EOP nutation correction in obliquity `δΔϵ`, in radians, 
at time `t` expressed in `TT` Julian centuries since J2000 for the IERS convention `m`.
"""
function eop_δΔϵ(::IERS2010A, t::Number)
    eop_check_init()
    return interpolate(IERS_EOP.nut2010.δΔϵ, t)
end

function eop_δΔϵ(::IERS2003A, t::Number)
    eop_check_init()
    return interpolate(IERS_EOP.nut2003.δΔϵ, t)
end

function eop_δΔϵ(::IERS1996, t::Number)
    eop_check_init()
    return interpolate(IERS_EOP.nut1996.δΔϵ, t)
end

eop_δΔϵ(::IERSModel, ::Number) = 0


"""
    eop_δX(m::IERSModel, t::Number)

Interpolate and retrieve the CIP `δX` correction, in radians, at time `t` expressed in 
`TT` Julian centuries since J2000 for the IERS convention `m`.
"""
function eop_δX(::IERS2010A, t::Number)
    eop_check_init()
    return interpolate(IERS_EOP.nut2010.δX, t)
end

function eop_δX(::IERS2003A, t::Number)
    eop_check_init()
    return interpolate(IERS_EOP.nut2003.δX, t)
end

function eop_δX(::IERS1996, t::Number)
    eop_check_init()
    return interpolate(IERS_EOP.nut1996.δX, t)
end

eop_δX(::IERSModel, ::Number) = 0


"""
    eop_δY(m::IERSModel, t::Number)

Interpolate and retrieve the CIP `δY` correction, in radians, at time `t` expressed in 
`TT` Julian centuries since J2000 for the IERS convention `m`.
"""
function eop_δY(::IERS2010A, t::Number)
    eop_check_init()
    return interpolate(IERS_EOP.nut2010.δY, t)
end

function eop_δY(::IERS2003A, t::Number)
    eop_check_init()
    return interpolate(IERS_EOP.nut2003.δY, t)
end

function eop_δY(::IERS1996, t::Number)
    eop_check_init()
    return interpolate(IERS_EOP.nut1996.δY, t)
end

eop_δY(::IERSModel, ::Number) = 0


"""
    eop_xp(m::IERSModel, t::Number)

Interpolate and retrieve the pole `xₚ` coordinate, in radians, at time `t` expressed in 
`TT` Julian centuries since J2000 for the IERS convention `m`.
"""
function eop_xp(::IERSModel, t::Number)
    eop_check_init()
    return interpolate(IERS_EOP.xp, t)
end

eop_xp(::CPND, ::Number) = 0


"""
    eop_yp(m::IERSModel, t::Number)

Interpolate and retrieve the pole `yₚ` coordinate, in radians, at time `t` expressed in 
`TT` Julian centuries since J2000 for the IERS convention `m`.
"""
function eop_yp(::IERSModel, t::Number)
    eop_check_init()
    return interpolate(IERS_EOP.yp, t)
end

eop_yp(::CPND, ::Number) = 0


"""
    offset_tt2ut1(t::Number)

Return the TT-to-UT1 offset, in seconds, at `t` expressed in TT seconds since J2000.
"""
function offset_tt2ut1(t)
    eop_check_init()
    return interpolate(IERS_EOP.ut1_tt, t/Tempo.CENTURY2SEC)
end


