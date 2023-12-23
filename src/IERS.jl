module IERS

    using Tempo

    using PrecompileTools
    using ReferenceFrameRotations
    using StaticArrays

    # Basic definitions
    include("conventions.jl")
    include("angles.jl")
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
    include("eop.jl")

    # Rotation functions 
    include("rotations.jl")

end