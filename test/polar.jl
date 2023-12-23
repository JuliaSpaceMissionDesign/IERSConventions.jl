
r2a = 180 / π * 3600
v2as = (x, y) -> acosd(max(-1, min(1, dot(x / norm(x), y / norm(y))))) * 3600


@testset "Polar Motion" verbose=true begin 

    for _ in 1:50

        m = rand((iers1996, iers2003a, iers2003b, iers2010a, iers2010b, CPNc, CPNd))

        tt_c = rand()
        tt_d = tt_c*Tempo.CENTURY2DAY

        v = rand(BigFloat, 3)

        # -- Testing TIO Locator 
        sp = IERS.tio_locator(iers2010a, tt_c) 
        spₑ = sp00(DJ2000, tt_d) 

        @test r2a*abs(sp-spₑ) ≤ 1e-6

        @test IERS.tio_locator(iers1996, tt_c) == 0
        @test IERS.tio_locator(CPNd, tt_c) == 0

        # -- Testing Polar Motion 
        xₚ, yₚ = rand(), rand()

        RPₐ = iers_polar_motion(iers2010a, xₚ, yₚ, tt_c)
        RPₑ = pom00(xₚ, yₚ, sp00(DJ2000, tt_d))

        @test v2as(RPₑ * v, RPₐ * v) ≤ 1e-8

        @test iers_polar_motion(CPNd, rand(), rand(), tt_c) == DCM(1I)

    end
end;