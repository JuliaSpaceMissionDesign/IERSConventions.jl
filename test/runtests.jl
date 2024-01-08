
using IERSConventions
using Test 

using DelimitedFiles
using ERFA 
using LinearAlgebra
using ReferenceFrameRotations
using RemoteFiles
using Tempo

# Download EOP files
EOP_DATA_FILE = @RemoteFile "https://datacenter.iers.org/data/csv/finals2000A.data.csv" dir = joinpath(
    @__DIR__, "assets"
);

download(EOP_DATA_FILE; verbose=true, force=false)

@info "Initialise EOP data"
let
    eopfile = joinpath(@__DIR__, "assets", "iau2000a")
    eop_generate_from_csv(iers2010a, path(EOP_DATA_FILE), eopfile)
    eop_load_data!(eopfile*".eop.dat", iers2010a)
end;

@testset "IERSConventions" verbose=true begin 
    
    include("fa.jl")
    include("bpn.jl")
    include("cio.jl")
    include("sidereal.jl")
    include("polar.jl")
    include("era.jl")

    # include("rotations.jl")

end;