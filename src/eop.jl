
# EOP Data structures 
# =========================================

# TODO: write me 


# EOP Retrival functions 
# =========================================

"""
    eop_δΔψ(m::IERSConventions, t::Number)

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

eop_δΔψ(::IERSConventions, ::Number) = 0


"""
    eop_δΔϵ(m::IERSConventions, t::Number)

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

eop_δΔϵ(::IERSConventions, ::Number) = 0


"""
    eop_δX(m::IERSConventions, t::Number)

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

eop_δX(::IERSConventions, ::Number) = 0


"""
    eop_δY(m::IERSConventions, t::Number)

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

eop_δY(::IERSConventions, ::Number) = 0


"""
    eop_xp(m::IERSConventions, t::Number)

Interpolate and retrieve the pole `xₚ` coordinate, in radians, at time `t` expressed in 
`TT` Julian centuries since J2000 for the IERS convention `m`.
"""
function eop_xp(m::IERSConventions, t::Number)
    # TODO: complete me
    return 0
end

eop_xp(::CPND, ::Number) = 0


"""
    eop_yp(m::IERSConventions, t::Number)

Interpolate and retrieve the pole `yₚ` coordinate, in radians, at time `t` expressed in 
`TT` Julian centuries since J2000 for the IERS convention `m`.
"""
function eop_yp(m::IERSConventions, t::Number)
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


# EOP Conversion functions 
# =========================================

# Function to convert nutation corrections to CIP corrections and viceversa
function δnut_to_δcip(m::IERSConventions, t::Number, δΔψ::Number, δΔϵ::Number)
    
    # Compute the precession angles 
    ϵ₀, ψₐ, _, χₐ = precession_angles_rot4(m, t)

    # Compute sine\cosine of mean obliquity
    se = sin(orient_obliquity(m, t))
    ce = cos(ϵ₀)

    c = ψₐ*ce - χₐ

    # Convert nutation corrections 
    δx = δΔψ*se + c*δΔϵ
    δy = δΔϵ - c*se*δΔψ

    return δx, δy

end

function δcip_to_δnut(m::IERSConventions, t::Number, δx::Number, δy::Number)

    # Compute the precession angles 
    ϵ₀, ψₐ, _, χₐ = precession_angles_rot4(m, t)

    # Compute sine\cosine of mean obliquity
    se = sin(orient_obliquity(m, t))
    ce = cos(ϵ₀)

    c = ψₐ*ce - χₐ
    d = 1 + c^2 

    # Convert CIP corrections 
    δΔψ = (δx - c*δy)/se/d
    δΔϵ = (δy + c*δx)/d

    return δΔψ, δΔϵ

end
