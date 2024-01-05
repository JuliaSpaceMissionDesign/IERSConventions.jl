
"""
    eop_δΔψ(m::IERSModel, t::Number)

Interpolate and retrieve the EOP nutation correction in longitude `δΔψ`, in radians, 
at time `t` expressed in `TT` Julian centuries since J2000 for the IERS convention `m`.
"""
function eop_δΔψ(m::IERSAModels, t::Number)
    # TODO: complete me
    return 0
end

function eop_δΔψ(m::IERS1996, t::Number)
    # TODO complete me 
    return 0
end

eop_δΔψ(::IERSModel, ::Number) = 0


"""
    eop_δΔϵ(m::IERSModel, t::Number)

Interpolate and retrieve the EOP nutation correction in obliquity `δΔϵ`, in radians, 
at time `t` expressed in `TT` Julian centuries since J2000 for the IERS convention `m`.
"""
function eop_δΔϵ(m::IERSAModels, t::Number)
    # TODO: complete me
    return 0
end

function eop_δΔϵ(m::IERS1996, t::Number)
    # TODO complete me 
    return 0
end

eop_δΔϵ(::IERSModel, ::Number) = 0


"""
    eop_δX(m::IERSModel, t::Number)

Interpolate and retrieve the CIP `δX` correction, in radians, at time `t` expressed in 
`TT` Julian centuries since J2000 for the IERS convention `m`.
"""
function eop_δX(m::IERSAModels, t::Number)
    # TODO: complete me
    return 0
end

function eop_δX(m::IERS1996, t::Number)
    # TODO complete me 
    return 0
end

eop_δX(::IERSModel, ::Number) = 0


"""
    eop_δY(m::IERSModel, t::Number)

Interpolate and retrieve the CIP `δY` correction, in radians, at time `t` expressed in 
`TT` Julian centuries since J2000 for the IERS convention `m`.
"""
function eop_δY(m::IERSAModels, t::Number)
    # TODO: complete me
    return 0
end

function eop_δY(m::IERS1996, t::Number)
    # TODO complete me 
    return 0
end

eop_δY(::IERSModel, ::Number) = 0


"""
    eop_xp(m::IERSModel, t::Number)

Interpolate and retrieve the pole `xₚ` coordinate, in radians, at time `t` expressed in 
`TT` Julian centuries since J2000 for the IERS convention `m`.
"""
function eop_xp(m::IERSModel, t::Number)
    # TODO: complete me
    return 0
end

eop_xp(::CPND, ::Number) = 0


"""
    eop_yp(m::IERSModel, t::Number)

Interpolate and retrieve the pole `yₚ` coordinate, in radians, at time `t` expressed in 
`TT` Julian centuries since J2000 for the IERS convention `m`.
"""
function eop_yp(m::IERSModel, t::Number)
    # TODO: complete me
    return 0
end

eop_yp(::CPND, ::Number) = 0


"""
    eop_Δut1(t::Number)

Interpolate and retrieve the `ΔUT1` value, in seconds, at time `t` expressed in `TT` Julian 
centuries since J2000, for the IERS convention `m`. 
"""
function eop_Δut1(t::Number)
    return 0 
end


function offset_tt2ut1(seconds)

    # Retrieve Leapseconds
    dAT = Tempo.offset_tai2utc(seconds - 32.184)

    # Retrieve ΔUT1 
    dUT1 = eop_Δut1(seconds/Tempo.CENTURY2SEC)

    return dUT1 + dAT - 32.184  

end


