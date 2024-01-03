
using IERSConventions
using Test 

using ERFA 
using LinearAlgebra
using ReferenceFrameRotations
using Tempo

@testset "IERSConventions" verbose=true begin 
    
    include("fa.jl")
    include("bpn.jl")
    include("cio.jl")
    include("sidereal.jl")
    include("polar.jl")
    include("era.jl")

    include("rotations.jl")

end;