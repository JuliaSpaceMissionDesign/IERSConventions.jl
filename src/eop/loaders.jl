
export eop_load_data!, eop_unload_data!

"""
    eop_load_data!(m::IERS2010, filename)

Initialise the Earth Orientation Parameters (EOP) from a dedicated JSMD `.eop.dat` file.
"""
function eop_load_data!(m::IERSModel, filename)

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
    IERS_EOP.nut1996 = NutCorrectionsInterpolator(IERS_EOP_DATA.cent_TT, IERS_EOP_DATA.nut1996)
    IERS_EOP.nut2003 = NutCorrectionsInterpolator(IERS_EOP_DATA.cent_TT, IERS_EOP_DATA.nut2003)
    IERS_EOP.nut2010 = NutCorrectionsInterpolator(IERS_EOP_DATA.cent_TT, IERS_EOP_DATA.nut2010)

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

    # Compute TT centuries
    tt_c = days_tt / Tempo.CENTURY2DAY

    # Set data time stamps
    IERS_EOP_DATA.filename = filename
    IERS_EOP_DATA.days_UTC = days_utc
    IERS_EOP_DATA.cent_TT = tt_c

    # Set TT-to-UT1 offset and LOD
    IERS_EOP_DATA.UT1_TT = ut1_tt
    IERS_EOP_DATA.LOD = lod

    # Set polar motion coordinates 
    IERS_EOP_DATA.xp = xp
    IERS_EOP_DATA.yp = yp

    # Transforma the CIP/nutation corrections from model `m` to the remaining ones.
    set_cip_nutation_corr!(m, tt_c, δX, δY, δΔψ, δΔϵ)
    nothing

end

function set_cip_nutation_corr!(m::IERSModel, tt_c, δX, δY, δΔψ, δΔϵ)

    n = length(tt_c)

    # Initialise CIP/Nutation corrections to the other models 
    δX_96, δY_96 = zeros(n), zeros(n)
    δX_03, δY_03 = zeros(n), zeros(n)
    δX_10, δY_10 = zeros(n), zeros(n)

    δΔψ_96, δΔϵ_96 = zeros(n), zeros(n)
    δΔψ_03, δΔϵ_03 = zeros(n), zeros(n)
    δΔψ_10, δΔϵ_10 = zeros(n), zeros(n)

    # radians to arcseconds factor
    k = 648000 / π

    @inbounds for j in eachindex(tt_c)

        # Compute CIP coordinates with the input model
        xin, yin = cip_xy(m, tt_c[j])

        # Retrieve corrections in radians
        δXj = δX[j] / k
        δYj = δY[j] / k

        δΔψj = δΔψ[j] / k
        δΔϵj = δΔϵ[j] / k

        # Converting CIP corrections to all conventions
        δX_96[j], δY_96[j] = convert_cip_corr(m, iers1996, tt_c[j], xin, yin, δXj, δYj)
        δX_03[j], δY_03[j] = convert_cip_corr(m, iers2003a, tt_c[j], xin, yin, δXj, δYj)
        δX_10[j], δY_10[j] = convert_cip_corr(m, iers2010a, tt_c[j], xin, yin, δXj, δYj)

        # Convert nutation corrections to all conventions
        δΔψ_96[j], δΔϵ_96[j] = convert_nut_corr(m, iers1996, tt_c[j], δX_96[j], δY_96[j], δΔψj, δΔϵj)
        δΔψ_03[j], δΔϵ_03[j] = convert_nut_corr(m, iers2003a, tt_c[j], δX_03[j], δY_03[j], δΔψj, δΔϵj)
        δΔψ_10[j], δΔϵ_10[j] = convert_nut_corr(m, iers2010a, tt_c[j], δX_10[j], δY_10[j], δΔψj, δΔϵj)

    end

    IERS_EOP_DATA.nut1996 = NutationCorrections(δX_96 * k, δY_96 * k, δΔψ_96 * k, δΔϵ_96 * k)
    IERS_EOP_DATA.nut2003 = NutationCorrections(δX_03 * k, δY_03 * k, δΔψ_03 * k, δΔϵ_03 * k)
    IERS_EOP_DATA.nut2010 = NutationCorrections(δX_10 * k, δY_10 * k, δΔψ_10 * k, δΔϵ_10 * k)

    nothing

end

function convert_cip_corr(::IERSModel, mout::IERSModel, tt_c, xin, yin, δX, δY)

    xout, yout = cip_xy(mout, tt_c)

    δX_out = xin - xout + δX
    δY_out = yin - yout + δY

    return δX_out, δY_out

end

function convert_cip_corr(::M, ::M, tt_c, x, y, δX, δY) where {M<:IERSModel}
    return δX, δY
end


function convert_nut_corr(::IERSModel, ::IERSModel, tt_c, δX, δY, δΔψ, δΔϵ)
    return δcip_to_δnut(iers2010a, tt_c, δX, δY)
end

function convert_nut_corr(::M, ::M, tt_c, δX, δY, δΔψ, δΔϵ) where {M<:IERSModel}
    return δΔψ, δΔϵ
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
                "the data with \'eop_generate_from_csv\' or \'eop_generate_from_txt\' and retry."
            )
        )
    end

    # Load the file
    data = readdlm(filename; header=false)

    # Retrieve UT1-UTC 
    Δut1 = @view(data[:, 4])

    # Convert UTC days to TT seconds
    days_utc = @view(data[:, 1])
    days_tai = map(t -> Tempo.utc2tai(Tempo.DJ2000, t)[2], days_utc)
    days_tt = days_tai .+ Tempo.OFFSET_TAI_TT / Tempo.DAY2SEC

    # Compute UT1-TT, in seconds
    sec_ut1 = days_utc * Tempo.DAY2SEC + Δut1
    ut1_tt = sec_ut1 - days_tt * Tempo.DAY2SEC

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
