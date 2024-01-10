
export eop_filename

# EOP Data 
# ============================

"""
    NutationCorrections

Container to hold the nutation corrections and celestial pole offsets associated to a  
given IERS model (e.g., 2000A)

### Fields    
- `δX, δY`: Celestial pole offsets referred to the model IAU2000A (rad).
- `δΔψ, δΔϵ`: Nutation corrections in longitude and obliquity (rad)

### See also 
See also [`EOPData`](@ref).
"""
struct NutationCorrections{T}
    δX::Vector{T}
    δY::Vector{T}
    δΔψ::Vector{T}
    δΔϵ::Vector{T}
end

function NutationCorrections(::Type{T}=Float64) where T 
    return NutationCorrections(T[], T[], T[], T[])
end


"""
    EOPData{T}

EOP Data container. Data is parameterised by time expressed in Terrestrial Time (TT) days. 

### Fields 
- `filename` : File where the EOP data are stored. 
- `days_UTC`: UTC Julian days since J2000.
- `cent_TT` : TT Julian centuries since J2000.
- `UT1_TT`: UT1 minus TT offset, in seconds.
- `LOD`: Length of day offset (s).
- `xp, yp`: Polar motion with respect to the crust (rad).
- `nut1996, nut2003, nut2006`: nutation and CIP corrections.  
"""
mutable struct EOPData{T}

    filename::String 

    days_UTC::Vector{T}
    cent_TT::Vector{T}

    UT1_TT::Vector{T}
    LOD::Vector{T}

    xp::Vector{T}
    yp::Vector{T}

    nut1996::NutationCorrections{T}
    nut2003::NutationCorrections{T}
    nut2010::NutationCorrections{T}

end

function EOPData(::Type{T}=Float64) where T 
    return EOPData(
        "", T[], [-Inf, -Inf], T[], T[], T[], T[], 
        NutationCorrections(T), NutationCorrections(T), NutationCorrections(T)
    )
end 

function Base.show(io::IO, eop::EOPData)

    if isempty(eop.days_UTC)
        println(io, "EOPData()")
    else 
        println(
            io, 
            "EOPData(filename=\"$(eop.filename)\", from: $(eop.days_UTC[1]) (UTC) to "* 
            "$(eop.days_UTC[end]) (UTC))"
        )
    end 

end

"""
    set_nutation_corr!(eop::EOPData, m::IERSModel, nc::NutationCorrections)

Set the nutation and cip corrections associated to the IERS model `m` to `nc`.
"""
set_nutation_corr!(eop::EOPData, ::IERS1996, nc::NutationCorrections) = eop.nut1996 = nc
set_nutation_corr!(eop::EOPData, ::IERS2003, nc::NutationCorrections) = eop.nut2003 = nc
set_nutation_corr!(eop::EOPData, ::IERS2010, nc::NutationCorrections) = eop.nut2010 = nc

"""
    eop_filename() 

Get the loaded Earth Orientation Parameters (EOP) filename.
"""
function eop_filename()
    
    if isempty(IERS_EOP_DATA.filename) 
        throw(ErrorException("Unable to retrieve filename, no EOP data has been loaded."))
    end

    return IERS_EOP_DATA.filename
end


# Interpolators 
# ============================

"""
    NutCorrectionsInterpolator

Container to store the interpolators for the nutation corrections and celestial pole offsets 
associated to a given IERS model. 

### See also 
See also [`EOPInterpolator`](@ref)
"""
struct NutCorrectionsInterpolator{T <: AbstractInterpolationMethod}
    δX::T
    δY::T
    δΔψ::T
    δΔϵ::T
end

function NutCorrectionsInterpolator(::Type{T}) where T
    NutCorrectionsInterpolator(
        _akima_init(T), _akima_init(T), _akima_init(T), _akima_init(T)
    )
end 

function NutCorrectionsInterpolator(t, nc::NutationCorrections) 
    NutCorrectionsInterpolator(
        InterpAkima(t, nc.δX), 
        InterpAkima(t, nc.δY), 
        InterpAkima(t, nc.δΔψ), 
        InterpAkima(t, nc.δΔϵ)
    )
end

"""
    EOPInterpolator

Container to store the interpolators for the loaded EOP data. 
"""
mutable struct EOPInterpolator{T <: AbstractInterpolationMethod}

    init::Bool 

    xp::T 
    yp::T 

    ut1_tt::T 
    lod::T
    
    nut1996::NutCorrectionsInterpolator{T}
    nut2003::NutCorrectionsInterpolator{T}
    nut2010::NutCorrectionsInterpolator{T}

end

function EOPInterpolator(::Type{T}=Float64) where T 
    EOPInterpolator(
        false, 
        _akima_init(T), _akima_init(T), _akima_init(T), _akima_init(T),
        NutCorrectionsInterpolator(T), 
        NutCorrectionsInterpolator(T), 
        NutCorrectionsInterpolator(T)
    )
end

function Base.show(io::IO, i::EOPInterpolator)
    println(io, "EOPInterpolator(init=$(i.init))")
end

# Checks whether EOP data has been initialised successfully
function eop_check_init()
    if !IERS_EOP.init
        throw(
            ErrorException(
                "EOP not initialized. Please run 'eop_load_data!' before using this function."
            )
        )
    end
    nothing
end

# Function to initialise a dummy Akima interpolator 
_akima_init(::Type{T}) where T = InterpAkima([1, 2, 3, 4, 5, 6], zeros(T, 6))


# IERS Constants  
# ============================

"""
    IERS_EOP_DATA 

Earth Orientation Parameters Data. 

### See also 
See also [`EOPData`](@ref).
"""
const IERS_EOP_DATA = EOPData();

"""
    IERS_EOP 

Earth Orientation Parameters (EOP) interpolators. 

### See also 
See also [`eop_load_data!`](@ref) and [`EOPInterpolator`](@ref).
"""
const IERS_EOP = EOPInterpolator();