@testset "Earth rotation rate" begin 
   @test IERSConventions.iers_earth_rot_rate(0.0) ≈ IERSConventions.iers_earth_rot_rate()
end