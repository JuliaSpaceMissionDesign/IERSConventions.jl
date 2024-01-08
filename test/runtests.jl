
using IERSConventions
using Test 

using DelimitedFiles
using ERFA 
using LazyArtifacts
using LinearAlgebra
using ReferenceFrameRotations
using RemoteFiles
using Tempo

# Download EOP files
EOP_DATA_FILE = @RemoteFile "https://datacenter.iers.org/data/csv/finals2000A.data.csv" dir = joinpath(
    @__DIR__, "assets"
);

download(EOP_DATA_FILE; verbose=true, force=false)

@testset "Download all artifacts" begin
    @info artifact"testdata"
    @info "All artifacts downloaded"
end;

test_dir = artifact"testdata"

begin 
    @info "Initialise EOP data"
    eopfile = joinpath(test_dir, "eopc04_20.1962-now.txt")
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

end;