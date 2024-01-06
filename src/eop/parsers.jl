
export eop_parse_csv, eop_parse_txt

"""
    eop_parse_csv(m::IERSModel, inputfile, outputfile)

Parse CSV files containing IERS EOP data and extracts the relevant information to a dedicated 
JSMD `.eop.dat` file. Supported formats are the EOP C04 series and the Rapid Data 
prediction (finals).

!!! note 
    The `outputfile` name should not include the file extension, which is automatically 
    added by this function.  

!!! note 
    Depending on the type of file, either the CIP or the nutation corrections may be present. 
    This function applies a conversion algorithm to automatically retrieve the missing data.

!!! warning 
    The rapid data prediction files store the LOD, CIP and nutation corrections in 
    milliseconds and milliarcseconds, respectively, whereas the EOP C04 files do not. This 
    routine automatically detects the correct unit of measure by analysing the format of the 
    input filename. Thus, the filename should be kept equal to the one released by the IERS. 

### References 
# TODO: add websites

### See also 
See also [`eop_parse_txt`](@ref). 
"""
function eop_parse_csv(m::IERSModel, inputfile, outputfile)

    # Check whether the file is a Combined Series (C04) or the Rapid Data Series.
    isfinal = startswith(splitdir(inputfile)[2], "finals")

    # Multiplicative factor to bring all quantities to seconds/arcseconds, since in the 
    # rapid data series the LOD is given in milliseconds, and the dX, dY, dPsi and dEps 
    # corrections are in milliarcseconds. 
    fct = isfinal ? 1e-3 : 1.0 
    
    # Read the data in the file 
    data, labels = readdlm(inputfile, ';'; header=true)   
    headers = collect(String.(labels[:]))

    # Instead of relying on fixed column indexes, we retrieve the column ID associated 
    # to a specific quantity by checking the header strings. If the required string is 
    # unavailable, then the CSV file is not supported! 

    cols = ["MJD", "LOD", "UT1-UTC", "x_pole", "y_pole", "dX", "dY", "dPsi", "dEpsilon"]

    idxs = zeros(Int, length(cols))
    for (j, label) in enumerate(cols) 
        # Check whether a valid index has been found! 
        idx  = findfirst(headers .== label)
        if !isnothing(idx)
            idxs[j] = idx
        else 
            throw(
                ErrorException(
                    "Unsupported or invalid EOP file. Could not find the '$(label)' column."
                )
            )
        end
    end

    # Retrieve the last row in which there is available data for ΔUT1, for predictions, 
    # the other parameters usually end before. 
    last_row = findlast(!isempty, @view(data[:, idxs[3]]))

    mjd  = convert(Vector{Float64}, @view(data[1:last_row, idxs[1]]))
    Δut1 = convert(Vector{Float64}, @view(data[1:last_row, idxs[3]]))

    # Convert UTC days to TT centuries (for the corrections conversion)
    days_utc = mjd .- Tempo.DMJD 
    days_tai = map(t->Tempo.utc2tai(Tempo.DJ2000, t)[2], days_utc)
    cent_tt  = (days_tai .+ Tempo.OFFSET_TAI_TT/Tempo.DAY2SEC)/Tempo.CENTURY2DAY
    
    # Retrieve pole coordinates 
    xp = convert(Vector{Float64}, @view(data[1:last_row, idxs[4]]))
    yp = convert(Vector{Float64}, @view(data[1:last_row, idxs[5]]))

    # Retrieve columns with optional missing data 
    raw_lod = @view(data[1:last_row, idxs[2]])
    raw_δX  = @view(data[1:last_row, idxs[6]])
    raw_δY  = @view(data[1:last_row, idxs[7]])
    raw_δΔψ = @view(data[1:last_row, idxs[8]])
    raw_δΔϵ = @view(data[1:last_row, idxs[9]])

    # Fill missing LOD elements
    lod = fct*fill_eop_data(raw_lod)

    k = π/648000
    
    # Retrieve the latest rows with valid data  
    lrow_nut = findlast(!isempty, raw_δΔψ)
    if isnothing(lrow_nut)
        # Need to convert CIP into nutation corrections 
        δX = fct*fill_eop_data(raw_δX)
        δY = fct*fill_eop_data(raw_δY)
    
        corr = map((t, dx, dy)->δcip_to_δnut(m, t, dx*k, dy*k), cent_tt, δX, δY)
        δΔψ, δΔϵ = map(x->x[1]/k, corr), map(x->x[2]/k, corr)
    
    else 
        # Need to convert nutation into CIP corrections
        δΔψ = fct*fill_eop_data(raw_δΔψ)
        δΔϵ = fct*fill_eop_data(raw_δΔϵ)
    
        corr = map((t, dp, de)->δnut_to_δcip(m, t, dp*k, de*k), cent_tt, δΔψ, δΔϵ)
        δΔψ, δΔϵ = map(x->x[1]/k, corr), map(x->x[2]/k, corr)
    end
    
    # Write the EOP data to the desired file
    eop_write_data(
        hcat(
            days_utc, xp, yp, Δut1, 
            round.(lod; digits=7), round.(δX; digits=7), round.(δY; digits=7), 
            round.(δΔψ; digits=7), round.(δΔϵ; digits=7)
        ), 
        outputfile
    )

end


"""
    eop_parse_txt(m::IERSModel, inputfile, outputfile)

Parse TXT files containing IERS EOP C04 data and extracts the relevant information to a 
dedicated JSMD `.eop.dat` file. 

!!! note 
    The `outputfile` name should not include the file extension, which is automatically 
    added by this function.  

!!! note 
    Depending on the type of file, either the CIP or the nutation corrections may be present. 
    This function applies a conversion algorithm to automatically retrieve the missing data.

!!! warning 
    This routine recognises the file structure and column ordering by analysing the 
    input file name, which should be left equal to the one retrieved from the IERS website.

### See also 
See also [`eop_parse_csv`](@ref). 
"""
function eop_parse_txt(m::IERSModel, inputfile, outputfile)

    # Check that the starting file pattern matches the EOPC04 nomenclature
    filename = splitdir(inputfile)[2]
    if !startswith(filename, "eopc04")
        throw(ArgumentError("Unsupported EOP C04 filename!"))
    end

    # Retrieve ITRF version
    itrfv = parse(Int, filename[8:9])
    if itrfv == 14 

        hascip = occursin("IAU2000", filename)
        idxs = [4, 8, 7, 5, 6, 9, 10]

        # In this version, the EOP data starts at the 15th row
        data = readdlm(inputfile; header=false, skipstart=14)

    elseif itrfv == 20 
        
        hascip = true
        idxs = [5, 13, 8, 6, 7, 9, 10]
        
        # Parse the data file. In this version, the header is completely commented
        data = readdlm(inputfile; header=false, comments=true)

    else 
        throw(ArgumentError("Unsupported ITRF version."))
    end

    # In this case we don't have issues of missing data. 
    # Thus we immediately convert everything. 
    mjd = @view(data[:, idxs[1]])
    lod = @view(data[:, idxs[2]])
    Δut1 = @view(data[:, idxs[3]])

    # Convert UTC days to TT centuries (for the correction conversion)
    days_utc = mjd .- Tempo.DMJD 
    days_tai = map(t->Tempo.utc2tai(Tempo.DJ2000, t)[2], days_utc)
    cent_tt  = (days_tai .+ Tempo.OFFSET_TAI_TT/Tempo.DAY2SEC)/Tempo.CENTURY2DAY

    # Retrieve the pole coordinates 
    xp = @view(data[:, idxs[4]])
    yp = @view(data[:, idxs[5]])

    k = π/648000
        
    # Retrieve the CIP and nutation corrections
    if hascip
        δX = @view(data[:, idxs[6]])
        δY = @view(data[:, idxs[7]])

        # Convert CIP into nutation corrections 
        corr = map((t, dx, dy)->δcip_to_δnut(m, t, dx*k, dy*k), cent_tt, δX, δY)
        δΔψ, δΔϵ = map(x->x[1]/k, corr), map(x->x[2]/k, corr)

    else 
        δΔψ = @view(data[:, idxs[6]])
        δΔϵ = @view(data[:, idxs[7]])

        corr = map((t, dp, de)->δnut_to_δcip(m, t, dp*k, de*k), cent_tt, δΔψ, δΔϵ)
        δΔψ, δΔϵ = map(x->x[1]/k, corr), map(x->x[2]/k, corr)
    end

    # Write the EOP data to the desired file
    eop_write_data(
        hcat(
            days_utc, xp, yp, Δut1, 
            round.(lod; digits=7), round.(δX; digits=7), round.(δY; digits=7), 
            round.(δΔψ; digits=7), round.(δΔϵ; digits=7)
        ), 
        outputfile
    )

end


"""
    eop_write_data(data, output_filename)

Write the EOP data stored in the matrix `data` to a dedicated JSMD `.eop-dat` file.  

!!! note 
    The `output_filename` should not include the file extension, which is automatically 
    added by this function.  

### See also 
See also [`eop_parse_csv`](@ref) and [`eop_parse_txt`](@ref). 
"""
function eop_write_data(data, output_filename)
    writedlm(output_filename * ".eop.dat", data)
    @info "IERS EOP file converted to '$(output_filename).eop.dat'."
    nothing
end


function fill_eop_data(raw_data)

    # Fill all the raws with empty data with the latest valid EOP element and 
    # return the output as a Float64 vector.

    data = zeros(size(raw_data))

    # Find index of the latest raw with valid data 
    lrow = findlast(!isempty, raw_data)
    isnothing(lrow) && return data 

    # Fill the missing elements with the latest available data
    data[1:lrow] = raw_data[1:lrow]
    data[lrow+1:end] .= raw_data[lrow]
    
    return data 

end


# EOP Conversion functions 
# =========================================

"""
    δnut_to_δcip(m::IERSModel, t::Number, δΔψ::Number, δΔϵ::Number)

Convert nutation corrections in longitude `δΔψ` and obliquity `δΔϵ` to CIP corrections 
`δX`, `δY` at time `t` expressed in TT Julian centuries since J2000 for model `m`. All 
input and output EOP are given in radians. 

### References 
- IERS Technical Note No. [36](https://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn36.html) 
- Paris IERS observatory FTP server (UAI2000 package).

### See also 
See also [`δcip_to_δnut`](@ref).
"""
function δnut_to_δcip(m::IERSModel, t::Number, δΔψ::Number, δΔϵ::Number)
    
    # Compute the precession angles 
    ϵ₀, ψₐ, _, χₐ = precession_angles_rot4(m, t)

    # Compute sine\cosine of mean obliquity
    se = sin(iers_obliquity(m, t))
    ce = cos(ϵ₀)

    c = ψₐ*ce - χₐ

    # Convert nutation corrections 
    δx = δΔψ*se + c*δΔϵ
    δy = δΔϵ - c*se*δΔψ

    return δx, δy

end


"""
    δcip_to_δnut(m::IERSModel, t::Number, δX::Number, δY::Number)

Convert CIP corrections `δX`, `δY` to nutation corrections in longitude `δΔψ` and obliquity 
`δΔϵ` to  at time `t` expressed in TT Julian centuries since J2000 for model `m`. All 
input and output EOP are given in radians. 

### References 
- IERS Technical Note No. [36](https://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn36.html) 
- Paris IERS observatory FTP server (UAI2000 package).

### See also 
See also [`δnut_to_δcip`](@ref).
"""
function δcip_to_δnut(m::IERSModel, t::Number, δx::Number, δy::Number)

    # Compute the precession angles 
    ϵ₀, ψₐ, _, χₐ = precession_angles_rot4(m, t)

    # Compute sine\cosine of mean obliquity
    se = sin(iers_obliquity(m, t))
    ce = cos(ϵ₀)

    c = ψₐ*ce - χₐ
    d = 1 + c^2 

    # Convert CIP corrections 
    δΔψ = (δx - c*δy)/se/d
    δΔϵ = (δy + c*δx)/d

    return δΔψ, δΔϵ

end
