
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

test_dir = artifact"testdata"

begin 
    @info "Initialise EOP data"
    eopfile = joinpath(test_dir, "eop", "eopc04_20.1962-now.txt")
    eop_generate_from_txt(iers2010a, eopfile, joinpath(@__DIR__, "assets", "eopc04"))
    eop_load_data!(joinpath(@__DIR__, "assets", "eopc04.eop.dat"), iers2010a)
end

@testset "IERSConventions" verbose=true begin 
    
    include("fa.jl")
    include("bpn.jl")
    include("cio.jl")
    include("sidereal.jl")
    include("polar.jl")
    include("era.jl")

    include("rotations.jl")

    include("parsers.jl")

end;