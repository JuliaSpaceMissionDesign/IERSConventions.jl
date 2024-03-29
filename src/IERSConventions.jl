module IERSConventions

    using DelimitedFiles
    using PrecompileTools
    using ReferenceFrameRotations
    using StaticArrays

    using Tempo
    using JSMDInterfaces.Math: AbstractInterpolationMethod, interpolate
    using JSMDUtils.Math: InterpAkima, arcsec2rad, skew

    # Basic definitions
    include("models.jl")
    include("poisson.jl")

    # Fundamental arguments
    include("delaunay.jl")
    include("planetary.jl")

    # Bias-Precession-Nutation
    include("obliquity.jl")

    include("precession.jl")
    include("nutation.jl")
    include("bpn.jl")
    include("cip.jl")
    include("fw.jl")

    # Time rotations
    include("era.jl")
    include("sidereal.jl")

    # Polar motion
    include("polar.jl")

    # EOP functions 
    include("eop/data.jl")
    include("eop/parsers.jl")
    include("eop/loaders.jl")
    include("eop/retrieval.jl")

    # Rotation functions 
    include("rate.jl")
    include("rotations.jl")

    # Package precompilation
    include("precompile.jl")

    function __init__()
        if !Tempo.has_timescale(TIMESCALES, Tempo.timescale_id(UT1))
            Tempo.add_timescale!(TIMESCALES, UT1, offset_tt2ut1; parent=TT)
        end
    end

end