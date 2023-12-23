
""" 
    PoissonSeries 

This container stores the non-polynomial part of a Poisson series, i.e., the 
sine and cosine coefficients and the Delaunay's and Planetary arguments of the nutation 
theory multipliers. 
"""
struct PoissonSeries{T}

    # sine and cosine coefficients
    s::T 
    c::T 

    dm::SVector{5, Int} # Delaunay arguments multipliers 
    pm::SVector{9, Int} # Planetary arguments multipliers

end


# POISSON SERIES builder 
# ===========================

# This function is used to automatically generated optimised functions that 
# compute the series associated to the IAU constants.
function build_series(
    fname::Symbol, iers_convention::Symbol, trig::AbstractVector, poly=nothing; 
    enable_pargs=true, unit_factor=arcsec2rad(1e-6)
)

    # Maps for the Delaunay's and Planetary arguments 
    DA_MAPS = SVector(:M, :S, :F, :D, :Ω)
    PA_MAPS = SVector(:λ_Me, :λ_Ve, :λ_Ea, :λ_Ma, :λ_Ju, :λ_Sa, :λ_Ur, :λ_Ne, :pₐ) 

    FA_MAPS = Vector(DA_MAPS)
    enable_pargs && append!(FA_MAPS, PA_MAPS)

    # Copy the expressions since it will be manipulated 
    ctrig = deepcopy(trig)

    # Retrieve number of series, i.e., the number of returned outputs
    nout = length(ctrig)
    if isnothing(poly)
        poly = [zeros(length(ctrig[i])) for i = 1:nout]
    elseif length(poly) != nout
        throw(ArgumentError("`trig` and `poly` must have the same length."))
    end

    # The elements of cpoly must always be >= ctrig
    fbody = Expr(:block)

    # Retrieve unique sets of coefficients multiplying the Fundamental Arguments 
    if enable_pargs
        fa_list = map(
            x->vcat(x.dm, x.pm), unique(x->vcat(x.dm, x.pm), vcat(vcat(ctrig...)...))
        )

    else
        fa_list = map(
            x->vcat(x.dm), unique(x->vcat(x.dm), vcat(vcat(ctrig...)...))
        )
    end

    # Pre-computes all the Fundamental Arguments multiplications
    fad = Dict{Symbol, Vector{Int}}()
    for (m, fm) in enumerate(FA_MAPS)

        # Retrieve the unique coefficients multiplying the m-th fundamental argument.
        fad[fm] = filter(x-> x != 0, unique(map(x->x[m], fa_list)))

        # Skip the remaining if the argument is never used.
        len = length(fad[fm])
        len == 0 && continue

        # Register the expressions for that factor
        name = Expr(:(.), fm in DA_MAPS ? :da : :pa, QuoteNode(fm))
        mult::Vector{Expr} = [
            Expr(:(=), Symbol("$fm$k"), Expr(:call, :*, fad[fm][k], name)) for k = 1:len
        ]

        append!(fbody.args, mult)

    end

    # Final @evalpoly call expression (depends on number of outputs)
    # Automatically includes transformation from arcseconds to radians
    epoly = Expr(:call, :*, Expr(:macrocall, Symbol("@evalpoly"), :, :t), unit_factor)
    pcall = [Expr(:(=), Symbol("x$j"), deepcopy(epoly)) for j = 1:nout]

    # Assings all the polynomial contributions
    for j in 1:nout
        for i = 1:max(length(ctrig[j]), length(poly[j]))
            pval = i > length(poly[j]) ? 0.0 : poly[j][i]
            push!(fbody.args, Expr(:(=), Symbol("x$(j)x$i"), pval))
            push!(pcall[j].args[2].args[2].args, Symbol("x$(j)x$i"))
        end
    end

    # For each unique ARGUMENT expression 
    for fa_set in reverse(fa_list) 

        # Stores the expressions for this FA set
        setₑ = Vector{Expr}()

        cₙ, sₙ = false, false
        for (j, series) in enumerate(ctrig) 
            for (i, block) in enumerate(series)

                # There should always be single one! 
                idx = findfirst(
                    x-> (enable_pargs ? vcat(x.dm, x.pm) : x.dm) == fa_set, block
                )

                isnothing(idx) && continue

                # Remove this element to make future searches faster
                row = popat!(ctrig[j][i], idx)

                if row.s != 0 || row.c != 0 
                    sₙ = sₙ ? sₙ : row.s != 0
                    cₙ = cₙ ? cₙ : row.c != 0

                    sumₑ = Expr(:call, :(+))
                    row.s != 0 && push!(sumₑ.args, Expr(:call, :*, row.s, :sarg))
                    row.c != 0 && push!(sumₑ.args, Expr(:call, :*, row.c, :carg))
                    
                    push!(setₑ, Expr(:(+=), Symbol("x$(j)x$i"), sumₑ))
                end

            end
        end

        # Computes the sin\cos for the ARGUMENT expression
        if cₙ || sₙ
            # Retrieve the list of non-null FA  
            arg_list = Symbol[]
            for (j, fa) in enumerate(fa_set)
                if fa != 0
                    push!(
                        arg_list, 
                        Symbol("$(FA_MAPS[j])$(argmin(abs.(fad[FA_MAPS[j]] .- fa)))")
                    )
                end
            end

            # Reassigns the ARG variable 
            if !isempty(arg_list)
                push!(fbody.args, Expr(:(=), :arg,  Expr(:call, :(+), arg_list...)))

                # Compute the sin(ARG) and cos(ARG)
                sₙ && push!(fbody.args, Expr(:(=), :sarg, Expr(:call, :sin, :arg)))
                cₙ && push!(fbody.args, Expr(:(=), :carg, Expr(:call, :cos, :arg)))

                push!(fbody.args, setₑ...)
            end
        end
    end

    for j in 1:nout
        push!(fbody.args, pcall[j])
    end

    # Adds return types
    if nout > 1
        rtupl = Expr(:tuple)
        for j in 1:nout
            push!(rtupl.args, Symbol("x$j"))
        end

        push!(fbody.args, Expr(:return, rtupl))
    end

    # Assembles the function
    fcn = _assemble_function(iers_convention, fname, fbody, enable_pargs)
    return eval(fcn)

end

function _assemble_function(convention, fname, fbody, pargs)

    # Creates function call Expression
    fcall = Expr(:call, fname)

    push!(fcall.args, Expr(:(::), convention))   # Adds model
    push!(fcall.args, Expr(:(::), :t, :Number)) # Adds time 

    push!(fcall.args, Expr(:(::), :da, :DelaunayArgs))  # Adds DA
    pargs && push!(fcall.args, Expr(:(::), :pa, :PlanetaryArgs)) # Adds PA

    # Assembles Function head and body
    return Expr(:function, fcall, fbody)

end



# IERS Table parsers 
# ===========================

# NOTE: the functions reported hereafter are currently intended only for internal use. They 
# are the ones used to generate the IERS constants files for the nutation and CIO series 
# from the original IERS tables. 


abstract type PoissonDataParser end 

# Number of Delaunay's arguments
n_dargs(::PoissonDataParser) = 5
n_dargs(::Type{<:PoissonDataParser}) = 5 

# Number of Planetary arguments  
n_pargs(::PoissonDataParser) = 9 
n_pargs(::Type{<:PoissonDataParser}) = 9

# Regex pattern for real numbers 
function rgx_real(::Type{<:PoissonDataParser})
    return "[-+]?(?:(?:[0-9]+(?:\\.[0-9]*)?)|(?:\\.[0-9]+))(?:[eE][-+]?[0-9]+)?"
end

# Regex pattern for integer numbers
function rgx_int(::Type{<:PoissonDataParser})
    return "[-+]?[0-9]+"
end

""" 
    CIODataParser <: PoissonDataParser 

Parser for CIO-like data, in which the orders of the different coefficients are 
separated by a specific line header in the form of " j = 0  Nb of terms = 1306". 

### Fields 
- ncols: number of columns 
- c_dargs: first column of the Delaunay's arguments 
- c_pargs: first column of the planetary arguments 
- c_sin: column with the sine terms 
- c_cos: column with the cosine terms 
- linePattern: regex expression to parse the relevant data lines  

"""
struct CIODataParser <: PoissonDataParser 
    ncols::Int 
 
    c_dargs::Int 
    c_pargs::Int 

    c_sin::Int
    c_cos::Int

    linePattern::Regex 
end

function CIODataParser(ncols::Int; c_sin::Int=-1, c_cos::Int=-1, kwargs...)
    create_poisson_parser(CIODataParser, ncols, c_sin, c_cos; kwargs...)
end

# Retrieve the first column of the Delaunay's and planetary arguments
fc_dargs(p::CIODataParser) = p.c_dargs
fc_pargs(p::CIODataParser) = p.c_pargs

# Retrieve Regex pattern for relevant data lines
pattern(p::CIODataParser) = p.linePattern

""" 
    NutationDataParser <: PoissonDataParser 

Parser for nutation-like data, in which the orders of the different coefficients are 
stored on separate columns, e.g., in and out of phase. 

### Fields 
- ncols: number of columns 
- c_dargs: first column of the Delaunay's arguments 
- c_pargs: first column of the planetary arguments 
- c_sin: column with the sine terms 
- c_cos: column with the cosine terms 
- linePattern: regex expression to parse the relevant data lines  

"""
struct NutationDataParser <: PoissonDataParser
    ncols::Int 

    c_dargs::Int 
    c_pargs::Int 

    c_sin::Vector{Int}
    c_cos::Vector{Int}
    linePattern::Regex 
end

function NutationDataParser(ncols::Int; c_sin=[-1], c_cos=[-1], kwargs...)
    create_poisson_parser(NutationDataParser, ncols, c_sin, c_cos; kwargs...)
end

# Retrieve the first column of the Delaunay's and planetary arguments
fc_dargs(p::NutationDataParser) = p.c_dargs
fc_pargs(p::NutationDataParser) = p.c_pargs

# Retrieve Regex pattern for relevant data lines
pattern(p::NutationDataParser) = p.linePattern

function create_poisson_parser(T::Type, ncols::Int, c_sin, c_cos; fc_del=-1, fc_pla=-1)
    
    # Assemble the line pattern reader 
    fields = ["\\S+" for _ in 1:ncols]

    # Retrieve number of Delaunay and planetary arguments 
    N_DEL = n_dargs(T)
    N_PLA = n_pargs(T)

    # Add the Delaunay and planetary arguments columns
    for (col, nargs) in zip((fc_del, fc_pla), (N_DEL, N_PLA))
        if col > 0 && col <= ncols + 1 - nargs 
            fields[col:col+nargs-1] .= rgx_int(T)
        end
    end
    
    # Add sine/cosine columns and update the fields entries 
    col_sin, col_cos = add_sincos!(fields, c_sin, c_cos)

    # Create the Regex expression 
    fields_str = "" 
    for (j, f) in enumerate(fields)
        fields_str *= "($f)"
        fields_str *= j < ncols ? "\\s+" : "\\s*\$"
    end

    T(ncols, fc_del, fc_pla, col_sin, col_cos, Regex(fields_str))
    
end

function add_sincos!(fields, c_sin::Int, c_cos::Int)
    # Retrieve regex for real numbers
    RGX_REAL = rgx_real(PoissonDataParser)

    # Add sine/cosine columns
    for c in (c_sin, c_cos)
        if c > 0 
            fields[c] = RGX_REAL
        end
    end

    return c_sin, c_cos
end

function add_sincos!(fields, c_sin::AbstractVector, c_cos::AbstractVector)
    # Retrieve regex for real numbers
    RGX_REAL = rgx_real(PoissonDataParser)

    # Add sine/cosine columns
    for cols in (c_sin, c_cos)
        for c in cols 
            if c > 0 
                fields[c] = RGX_REAL
            end
        end
    end

    # Retrieve the sizes of both vectors 
    ls, lc = length(c_sin), length(c_cos)

    # Number of expected degrees (usually 2)
    ndegs = max(ls, lc) 

    col_sin = zeros(Int, ndegs)
    col_cos = zeros(Int, ndegs)

    # Store sin\cos columns 
    col_sin[1:ls] .= c_sin 
    col_cos[1:lc] .= c_cos

    return col_sin, col_cos
end

function parse_arguments(n::Int, fcol::Int, match)

    # If no arguments have been specified, return an empty vector 
    fcol < 1 && return SVector{n}(zeros(Int, n))

    # Retrieve the arguments 
    args = match.captures[fcol:fcol + n - 1]
    return SVector{n}([parse(Int, x) for x in args])

end

function add_series_entry!(data, p::CIODataParser, m)

    # Retrieve the Delaunay's arguments if specified 
    dm = parse_arguments(n_dargs(p), fc_dargs(p), m)

    # Retrieve the planetary arguments if specified 
    pm = parse_arguments(n_pargs(p), fc_pargs(p), m)

    # Retrieve the sin/cos factors 
    cc = p.c_cos > 0 ? parse(Float64, m.captures[p.c_cos]) : 0.0
    sc = p.c_sin > 0 ? parse(Float64, m.captures[p.c_sin]) : 0.0 

    # Add the new series to the overall data vector
    if cc != 0 || sc != 0
        # No point in adding it if both entries are null
        push!(data[end], PoissonSeries(sc, cc, dm, pm))
    end

    return 1

end

function add_series_entry!(data, p::NutationDataParser, m)

    # Retrieve the Delaunay's arguments if specified 
    dm = parse_arguments(n_dargs(p), fc_dargs(p), m)

    # Retrieve the planetary arguments if specified 
    pm = parse_arguments(n_pargs(p), fc_pargs(p), m)

    # Retrieve the sin/cos factors 
    for j in eachindex(p.c_cos)

        cc = p.c_cos[j] > 0 ? parse(Float64, m.captures[p.c_cos[j]]) : 0.0
        sc = p.c_sin[j] > 0 ? parse(Float64, m.captures[p.c_sin[j]]) : 0.0 

        # Add the new series to the overall data vector
        if cc != 0 || sc != 0 
            # No point in adding it if both entries are null
            push!(data[j], PoissonSeries(sc, cc, dm, pm))
        end
    end
    nothing

end

"""
    parse_iers_constants(filename::String, parser::PoissonDataParser)

Retrieve the Poisson series within `filename` with the configuration of `parser`.
"""
function parse_iers_constants(filename::String, parser::CIODataParser)

    # Regex expression for the header (required only for parsing CIO series)
    reg_header = r"\s*j\s*=\s*([0-9]+)\s*[A-Za-z\s]+=\s([0-9]+)\s*$"

    # Initialise the data vector 
    data = []

    # Current entry, degree and total terms of this order 
    crt, deg, nterms = 0, 0, 0

    # Parse the input file
    open(filename, "r") do file 
        for line in readlines(file)
            if occursin(reg_header, line)
    
                # Check that all previous terms have been succesfully loaded 
                if nterms > 0 && crt != nterms 
                    throw(
                        ErrorException(
                            "Not all terms of degree $deg have been read:"*
                            " expected $nterms, parsed $crt."
                        )
                    )
                end

                m = match(reg_header, line)
                
                # Time argument exponent
                deg = parse(Int, m.captures[1])

                # Number of terms for this order
                nterms = parse(Int, m.captures[2])

                # Reset the terms counter 
                crt = 0
                
                # Create a new entry in the vector
                push!(data, PoissonSeries[])
    
            elseif occursin(pattern(parser), line)
                m = match(pattern(parser), line)
                crt += add_series_entry!(data, parser, m)
            end        
        end
    end

    return data 

end

function parse_iers_constants(filename::String, parser::NutationDataParser)

    # Initialise the data vector 
    data = [[] for _ in eachindex(parser.c_cos)]

    # Parse the input file
    open(filename, "r") do file 
        for line in readlines(file)
            if occursin(pattern(parser), line)
                m = match(pattern(parser), line)
                add_series_entry!(data, parser, m)
            end        
        end
    end

    return data 
end


"""
    generate_iers_file(filename::String, constant_name::Symbol, constant_data)

Generate or append the data in `constant_data` to `filename` as a variable named 
`constant_name`.
"""
function generate_iers_file(filename::String, constant_name::Symbol, constant_data)

    # Compute total number of terms 
    ntot = 0 
    for set in constant_data
        ntot += length(set)
    end

    # Write the parsed data into a file 
    open(filename, "a") do f 

        write(f,  "\nconst $constant_name = [\n")

        prv = 0
        for (j, set) in enumerate(constant_data)
    
            write(f, "\n\t# j = $(j-1)\n")
            write(f, "\t[")

            for (k, s) in enumerate(set)

                # Add a small comment every 10 terms
                if mod(k, 10) == 1
                    write(f, "\n\t\t# $(k+prv)-$(min(prv+k+9, prv+length(set)))\n")
                end

                # Add the series 
                write(f, "\t\tPoissonSeries($(s.s), $(s.c)")
                
                for v in (s.dm, s.pm)
                    write(f, ", SVector")
                    write(f, "$(repr(tuple(v...)))")
                end 

                write(f, "),\n")

            end

            write(f, "\t],\n")
            prv += length(set)

        end

        write(f, "\n]\n")
    end

end
