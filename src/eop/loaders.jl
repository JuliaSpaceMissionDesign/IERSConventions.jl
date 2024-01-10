
export eop_load_data!, eop_unload_data!

"""
    eop_load_data!(m::IERS2010, filename)

Initialise the Earth Orientation Parameters (EOP) from a dedicated JSMD `.eop.dat` file.

!!! note 
    Currently, only EOP data associated to the IAU2006/2000A precession-nutation model 
    is supported.
"""
function eop_load_data!(m::IERS2010, filename)

    # Set the EOP data 
    eop_set_data!(m, filename)

    # Initialise the interpolators 
    IERS_EOP.init = true 

    # Set polar motion
    IERS_EOP.xp = InterpAkima(IERS_EOP_DATA.cent_TT, IERS_EOP_DATA.xp)
    IERS_EOP.yp = InterpAkima(IERS_EOP_DATA.cent_TT, IERS_EOP_DATA.yp)

    IERS_EOP.lod = InterpAkima(IERS_EOP_DATA.cent_TT, IERS_EOP_DATA.LOD)
    IERS_EOP.ut1_tt = InterpAkima(IERS_EOP_DATA.cent_TT, IERS_EOP_DATA.UT1_TT)
    
    # Set nutation correction interpolators
    # TODO: ensure these are set! 
    # nci1996 = NutCorrectionsInterpolator(IERS_EOP_DATA.cent_TT, IERS_EOP_DATA.nut1996)
    # nci2003 = NutCorrectionsInterpolator(IERS_EOP_DATA.cent_TT, IERS_EOP_DATA.nut2003)
    nci2010 = NutCorrectionsInterpolator(IERS_EOP_DATA.cent_TT, IERS_EOP_DATA.nut2010)

    # IERS_EOP.nut1996 = nci1996
    # IERS_EOP.nut2003 = nci2003 
    IERS_EOP.nut2010 = nci2010

    @info "EOP initialised from file `$(filename)`."
    nothing 

end

""" 
    eop_unload_data!()

Unload all the EOP data that has been loaded during the session. 
"""
function eop_unload_data!()
    
    # Here we are not actually removing the data from the EOP data structures, 
    # but we still make sure that the parameters that are used to verify whether 
    # there is available EOP data are set to false. 

    IERS_EOP_DATA.filename = ""
    IERS_EOP.init = false

    @info "EOP data successfully unloaded."
    nothing 

end


"""
    eop_set_data!(m::IERSModel, filename)  

Set Earth Orientation Parameters (EOP) to be used for frames transformations from a JSMD 
`.eop.dat` file, given a reference model `m`.
"""
function eop_set_data!(m::IERSModel, filename)

    oldfile = IERS_EOP_DATA.filename
    days_utc, days_tt, xp, yp, ut1_tt, lod, δX, δY, δΔψ, δΔϵ = eop_read_data(filename)

    if (!isempty(oldfile))
        @warn "Existing EOP data from \'$(oldfile)\' will be overwritten by \'$(filename)\'."
    end

    # Set data time stamps
    IERS_EOP_DATA.filename = filename
    IERS_EOP_DATA.days_UTC = days_utc
    IERS_EOP_DATA.cent_TT  = days_tt/Tempo.CENTURY2DAY 

    # Set TT-to-UT1 offset and LOD
    IERS_EOP_DATA.UT1_TT = ut1_tt 
    IERS_EOP_DATA.LOD    = lod

    # Set polar motion coordinates 
    IERS_EOP_DATA.xp = xp 
    IERS_EOP_DATA.yp = yp

    # Set the nutation corrections 
    nc = NutationCorrections(δX, δY, δΔψ, δΔϵ)
    set_nutation_corr!(IERS_EOP_DATA, m, nc) 

    # TODO: convert these corrections to those for the other models!
    nothing

end


"""
    eop_read_data(filename)  

Read Earth Orientation Parameters (EOP) from a dedicated JSMD `.eop.dat` file. 

!!! note 
    JSMD's `.eop.dat` EOP files are meant to have a fixed structure to ease the 
    retrieval of the relevant EOP data. 
"""
function eop_read_data(filename)
    
    if !endswith(filename, "eop.dat")
        throw(
            ArgumentError(
                "EOP reader support only '.eop.dat' files! Please prepare " * 
                "the data with \'eop_generate_from_csv\' and retry."
            )
        )
    end

    # Load the file
    data = readdlm(filename; header=false)

    # Retrieve UT1-UTC 
    Δut1 = @view(data[:, 4])
    
    # Convert UTC days to TT seconds
    days_utc = @view(data[:, 1])
    days_tai = map(t->Tempo.utc2tai(Tempo.DJ2000, t)[2], days_utc)
    days_tt  = days_tai .+ Tempo.OFFSET_TAI_TT/Tempo.DAY2SEC

    # Compute UT1-TT, in seconds
    sec_ut1 = days_utc*Tempo.DAY2SEC + Δut1
    ut1_tt = sec_ut1 - days_tt*Tempo.DAY2SEC

    # Retrieve other quantities
    xp = @view(data[:, 2])
    yp = @view(data[:, 3])
    lod = @view(data[:, 5])

    δX = data[:, 6]
    δY = data[:, 7]
    δΔψ = data[:, 8]
    δΔϵ = data[:, 9]

    return days_utc, days_tt, xp, yp, ut1_tt, lod, δX, δY, δΔψ, δΔϵ

end