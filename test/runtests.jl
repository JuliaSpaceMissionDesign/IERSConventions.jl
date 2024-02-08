
using IERSConventions
using Test 

using DelimitedFiles
using ERFA 
using LazyArtifacts
using LinearAlgebra
using ReferenceFrameRotations
using RemoteFiles
using StaticArrays
using Tempo


@testset "Download all artifacts" begin
    @info artifact"testdata"
    @info "All artifacts downloaded"
end;

@testset "IERSConventions" verbose=true begin 

    include("eop.jl")
    include("rate.jl")

    include("fa.jl")
    include("bpn.jl")
    include("cio.jl")
    include("sidereal.jl")
    include("polar.jl")
    include("era.jl")

    include("rotations.jl")

    include("parsers.jl")

end;