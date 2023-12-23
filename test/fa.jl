
r2a = 180 / π * 3600

@testset "Delaunay Arguments" verbose=true begin

    r2a = 180 / π * 3600

    da = [
        IERS.DelaunayArgs(iers1996, 0.06567),
        IERS.DelaunayArgs(iers1996, 1.0)
    ]

    @testset "1996" begin 
        # -- Testing Mean Anomaly of the Moon (< 1μas)
        @test r2a*abs(da[1].M - 2.66359306441736)                   ≤ 1e-6
        @test r2a*abs(da[2].M - mod(-4.56593937247721e-1, 2π))      ≤ 1e-6

        # -- Testing Mean Anomaly of the Sun (< 1μas)
        @test r2a*abs(da[1].S - mod(-2.76485707808279, 2π))         ≤ 1e-6
        @test r2a*abs(da[2].S - mod(-5.97269171806332e-2, 2π))      ≤ 1e-6

        # -- Testing Mean Argument of Latitude of the Moon (< 1μas)
        @test r2a*abs(da[1].F - 2.53331724178132)                   ≤ 1e-6
        @test r2a*abs(da[2].F - 3.05931379900095)                   ≤ 1e-6

        # -- Testing mean elongation of the moon from the sun (< 1μas)
        @test r2a*abs(da[1].D - 3.23611369830128e-1)                ≤ 1e-6
        @test r2a*abs(da[2].D - mod(-2.00782792050270, 2π))         ≤ 1e-6

        # -- Testing Mean Longitude of the Moon (< 1μas)
        @test r2a*abs(da[1].Ω - mod(-3.43864262297640e-02, 2π))     ≤ 1e-6
        @test r2a*abs(da[2].Ω - mod(-1.58644591849564e-01, 2π))     ≤ 1e-6
    end

    t = [rand(49)..., 1.0]
    @testset "2003/2010 A" begin 
        for i in eachindex(t)
            da = IERS.DelaunayArgs(rand((iers2003a, iers2010a)), t[i])

            # --- Delaunay Arguments (accurate to 1μas)
            @test r2a*(abs(da.M - fal03(t[i])))                 ≤ 1e-6
            @test r2a*(abs(da.S - falp03(t[i])))                ≤ 1e-6
            @test r2a*(abs(da.F - faf03(t[i])))                 ≤ 1e-6
            @test r2a*(abs(da.D - fad03(t[i])))                 ≤ 1e-6
            @test r2a*(abs(da.Ω - mod(faom03(t[i]), 2π)))       ≤ 1e-6

        end
    end

    # Testing the truncated expressions of the Delunary Arguments 
    dab = [
        IERS.DelaunayArgs(rand([iers2003b, iers2010b]), 0.06567),
        IERS.DelaunayArgs(rand([iers2003b, iers2010b]), 1.0)
    ]

    @testset "2003/2010 B" begin

        # -- Testing Mean Anomaly of the Moon (< 1μas)
        @test r2a*abs(dab[1].M - 2.663599945842266)             ≤ 1e-6
        @test r2a*abs(dab[2].M - 5.826449449627781)             ≤ 1e-6

        # -- Testing Mean Anomaly of the Sun (< 1μas)
        @test r2a*abs(dab[1].S - 3.518352372771510)             ≤ 1e-6
        @test r2a*abs(dab[2].S - 6.223484580361229)             ≤ 1e-6

        # -- Testing Mean Argument of Latitude of the Moon (< 1μas)
        @test r2a*abs(dab[1].F - 2.533320574432192)             ≤ 1e-6
        @test r2a*abs(dab[2].F - 3.059379762905646)             ≤ 1e-6

        # -- Testing mean elongation of the moon from the sun (< 1μas)
        @test r2a*abs(dab[1].D - 3.236085510636535e-1)          ≤ 1e-6
        @test r2a*abs(dab[2].D - 4.275387201216140)             ≤ 1e-6

        # -- Testing Mean Longitude of the Moon (< 1μas)
        @test r2a*abs(dab[1].Ω - -3.438601115926876e-2 - 2π)    ≤ 1e-6
        @test r2a*abs(dab[2].Ω - -1.586802211172697e-1 - 2π)    ≤ 1e-6

    end

end;

@testset "Planetary Arguments" verbose=true begin

    t = [rand(49)..., 1.0]

    for i in eachindex(t)
        pa = IERS.PlanetaryArgs(rand((iers2003a, iers2010a)), t[i])

        # --- Planetary Arguments (accurate to 1μas)
        @test r2a*(abs(pa.λ_Me - fame03(t[i]))) ≤ 1e-6
        @test r2a*(abs(pa.λ_Ve - fave03(t[i]))) ≤ 1e-6
        @test r2a*(abs(pa.λ_Ea - fae03(t[i]) )) ≤ 1e-6
        @test r2a*(abs(pa.λ_Ma - fama03(t[i]))) ≤ 1e-6
        @test r2a*(abs(pa.λ_Ju - faju03(t[i]))) ≤ 1e-6
        @test r2a*(abs(pa.λ_Sa - fasa03(t[i]))) ≤ 1e-6
        @test r2a*(abs(pa.λ_Ur - faur03(t[i]))) ≤ 1e-6
        @test r2a*(abs(pa.λ_Ne - fane03(t[i]))) ≤ 1e-6
        @test r2a*(abs(pa.pₐ   - fapa03(t[i]))) ≤ 1e-6

    end

    # Check empty constructor
    pa = IERS.PlanetaryArgs(iers2003b, rand())
    for j in 1:9
        @test getfield(pa, j) == 0
    end

    # Check error for 1996 models 
    for fcn in (
        IERS.pa_mercury, IERS.pa_uranus, IERS.pa_neptune
    )
    
        @test_throws ErrorException fcn(iers1996, 0.0)

    end

end;