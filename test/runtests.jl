
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
    eop_parse_csv(iers2010a, path(EOP_DATA_FILE), eopfile)
    eop_load_data!(eopfile*".eop.dat", iers2010a)
end;

@testset "IERSConventions" verbose=true begin 
    
    include("fa.jl")
    include("bpn.jl")
    include("cio.jl")
    include("sidereal.jl")
    include("polar.jl")
    include("era.jl")

    include("rotations.jl")

end;

tt_c = rand()/4
tt_d = tt_c*Tempo.CENTURY2DAY

ut1_d = tt_d + IERSConventions.offset_tt2ut1(tt_d*Tempo.DAY2SEC)/Tempo.DAY2SEC

# --- Testing 1996 model (< 10 μas)
# GMST 
gm  = iers_gmst(iers1996, tt_c)
gme = gmst82(DJ2000, ut1_d)
@test r2a*abs(gm-gme) ≤ 1e-5