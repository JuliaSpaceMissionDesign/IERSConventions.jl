
@testset "IERS" verbose = true begin

    atol, rtol = 1e-11, 1e-11

    # Radians to arcseconds  
    r2a = 180 / π * 3600

    # Function to compute the angle between 2 vectors in arcseconds 
    v2as = (x, y) -> acosd(max(-1, min(1, dot(x / norm(x), y / norm(y))))) * 3600

    # Testing Delaunay's Arguments 
    @testset "Delaunay Arguments" verbose=true begin

        da = [
            DelaunayArgs(iers1996, 0.06567),
            DelaunayArgs(iers1996, 1.0)
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
                da = DelaunayArgs(rand((iers2003a, iers2010a)), t[i])

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
            DelaunayArgs(rand([iers2003b, iers2010b]), 0.06567),
            DelaunayArgs(rand([iers2003b, iers2010b]), 1.0)
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

    end

    @testset "Planetary Arguments" verbose=true begin

        t = [rand(49)..., 1.0]

        for i in eachindex(t)
            pa = PlanetaryArgs(rand((iers2003a, iers2010a)), t[i])

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
        pa = PlanetaryArgs(iers2003b, rand())
        for j in 1:9
            @test getfield(pa, j) == 0
        end

        # Check error for 1996 models 
        for fcn in (
            pa_mercury, pa_uranus, pa_neptune
        )
        
            @test_throws ErrorException fcn(iers1996, 0.0)

        end

    end

    @testset "Obliquity" verbose = true begin
        atol, rtol = 1e-8, 1e-8
        
        for _ in 1:50

            tt_c = rand()/2
            tt_d = tt_c*Tempo.CENTURY2DAY

            # For the 2003 models we also add the precession-rate corrections 
            ϵ96 = obl80(DJ2000, tt_d)
            ϵ00 = ϵ96 + pr00(DJ2000, tt_d)[2]
            ϵ06 = obl06(DJ2000, tt_d)
            
            # Accurate up to 1nas
            @test r2a*(abs(orient_obliquity(iers1996, tt_c) - ϵ96))     ≤ 1e-9
            @test r2a*(abs(orient_obliquity(iers2003a, tt_c) - ϵ00))    ≤ 1e-9
            @test r2a*(abs(orient_obliquity(iers2010a, tt_c) - ϵ06))    ≤ 1e-9

        end
    end

    @testset "Nutation" verbose=true begin 
        @testset "Nutation Components" verbose = true begin
            for _ in 1:50
                tt_c = rand()/2 
                tt_d = tt_c*Tempo.CENTURY2DAY

                # -- Testing IERS 1996 model (< 0.1μas)
                Δψ, Δϵ = orient_nutation_comp(iers1996, tt_c)
                Δp, Δe = nut80(DJ2000, tt_d)

                @test r2a*abs(Δψ - Δp) ≤ 1e-7
                @test r2a*abs(Δϵ - Δe) ≤ 1e-7

                # -- Testing IERS 2000A model (< 0.1μas)
                Δψ, Δϵ = orient_nutation_comp(iers2003a, tt_c)
                Δp, Δe = nut00a(DJ2000, tt_d)

                @test r2a*abs(Δψ - Δp) ≤ 1e-7
                @test r2a*abs(Δϵ - Δe) ≤ 1e-7

                # -- Testing IERS 2000B model (< 1mas)
                Δψ, Δϵ = orient_nutation_comp(iers2003b, tt_c)
                Δp, Δe = nut00b(DJ2000, tt_d)

                @test r2a*abs(Δψ - Δp) ≤ 1e-7
                @test r2a*abs(Δϵ - Δe) ≤ 1e-7

                # -- Testing IERS 2010A model (< 0.1μas)
                Δψ, Δϵ = orient_nutation_comp(iers2010a, tt_c)
                Δp, Δe = nut06a(DJ2000, tt_d)

                @test r2a*abs(Δψ - Δp) ≤ 1e-7
                @test r2a*abs(Δϵ - Δe) ≤ 1e-7

                # -- Testing IERS 2010B model (< 2mas)
                Δψ, Δϵ = orient_nutation_comp(iers2010b, tt_c)
                @test r2a*abs(Δψ - Δp) ≤ 2e-3
                @test r2a*abs(Δϵ - Δe) ≤ 2e-3
                
                # -- Testing IERS CPNc model (< 30mas)
                Δψ, Δϵ = orient_nutation_comp(CPNc, tt_c)
                @test r2a*abs(Δψ - Δp) ≤ 35e-3
                @test r2a*abs(Δϵ - Δe) ≤ 35e-3

                # -- Testing IERS CPNd model (< 1as)
                Δψ, Δϵ = orient_nutation_comp(CPNd, tt_c)
                @test r2a*abs(Δψ - Δp) ≤ 1
                @test r2a*abs(Δϵ - Δe) ≤ 1
            end
        end

        @testset "Nutation Matrix" verbose=true begin 
            for _ in 1:50
                tt_c = rand() 
                tt_d = tt_c*Tempo.CENTURY2DAY

                v = randn(BigFloat, 3)
                v /= norm(v)

                # --- Testing IERS 1996 model (< 0.1μas)
                N1 = orient_nutation(iers1996, tt_c)
                Ne = nutm80(DJ2000, tt_d)
                @test v2as(N1*v, Ne*v) ≤ 1e-7

                # --- Testing IERS 2003 models

                Ne = num00a(DJ2000, tt_d)
                ve = Ne*v; ve /= norm(ve);

                # IERS 2000A model (< 0.1μas)
                N1 = orient_nutation(iers2003a, tt_c)
                @test v2as(N1*v, ve) ≤ 1e-7

                # IERS 2000B model (< 1mas)
                N1 = orient_nutation(iers2003b, tt_c)
                @test v2as(N1*v, ve) ≤ 2e-3

                # --- Testing IERS 2010 models

                Ne = num06a(DJ2000, tt_d)
                ve = Ne*v; ve /= norm(ve);

                # IERS 2010A model (< 0.1μas)
                N1 = orient_nutation(iers2010a, tt_c)
                @test v2as(N1*v, ve) ≤ 1e-7

                # IERS 2010B model (< 1mas)
                N1 = orient_nutation(iers2010b, tt_c)
                @test v2as(N1*v, ve) ≤ 2e-3

                # CPNc model (< 35mas)
                N1 = orient_nutation(CPNc, tt_c)
                @test v2as(N1*v, ve) ≤ 35e-3

                # CPNd model (< 1mas)
                N1 = orient_nutation(CPNd, tt_c)
                @test v2as(N1*v, ve) ≤ 1

            end

        end

    end

    @testset "Precession" verbose=true begin 

        for j in 1:50

            tt_c = rand()
            
            j == 1 && (tt_c = 1)
            tt_d = tt_c*Tempo.CENTURY2DAY

            v = randn(BigFloat, 3)
            v /= norm(v)

            # -- Testing 1996 models (< 1 μas)
            Pe = pmat76(DJ2000, tt_d)
            v1 = Pe*v; v1 /= norm(v1); 

            # Test traditional parameterization
            zₐ, θₐ, ζₐ = precession_angles_rot3(iers1996, tt_c)
            ζe, ze, θe = prec76(DJ2000, 0, DJ2000, tt_d)
            
            @test r2a*abs(zₐ - ze) ≤ 1e-6
            @test r2a*abs(θₐ - θe) ≤ 1e-6
            @test r2a*abs(ζₐ - ζe) ≤ 1e-6

            # Test Capitaine parameterization (since we do not have ERFA expressions 
            # that provide this angles, we check that the resulting precession matrix is 
            # the same as that obtained with the original parameterization)
            ϵ₀, ψₐ, ωₐ, χₐ = precession_angles_rot4(iers1996, tt_c)
            P = angle_to_dcm(χₐ, :Z)*angle_to_dcm(ϵ₀, -ψₐ, -ωₐ, :XZX)

            # This test is accurate up to 0.1 mas, which is anyway below the precision of 
            # the 1977 precession model. 
            @test v2as(v1, P*v) ≤ 1e-4

            # Precession matrix 
            Pi = orient_precession(iers1996, tt_c)
            @test v2as(v1, Pi*v) ≤ 1e-6

            # --- Testing IERS 2003 models (< 1 μas)
            m = rand((iers2003a, iers2003b))

            Pe = bp00(DJ2000, tt_d)[2]
            v1 = Pe*v; v1 /= norm(v1); 

            # Test traditional parameterization
            zₐ, θₐ, ζₐ = precession_angles_rot3(iers2003a, tt_c)
            P = angle_to_dcm(-ζₐ, θₐ, -zₐ, :ZYZ)
            @test v2as(v1, P*v) ≤ 1e-6

            # Precession matrix 
            # This one also tests the Capitaine parameterization. 
            Pi = orient_precession(m, tt_c)
            @test v2as(v1, Pi*v) ≤ 1e-6

            # --- Testing IERS 2010 models (< 1 μas)
            m = rand((iers2010a, iers2010b, CPNc, CPNd))

            # Retrieve all the angles 
            ϵe, ψe, ωe, χe, ze, ζe, θe = p06e(DJ2000, tt_d)[[1, 2, 3, 9, 10, 11, 12]]

            # Test equatorial precession angles 
            zₐ, θₐ, ζₐ = precession_angles_rot3(m, tt_c)

            @test r2a*abs(zₐ - ze) ≤ 1e-6
            @test r2a*abs(θₐ - θe) ≤ 1e-6
            @test r2a*abs(ζₐ - ζe) ≤ 1e-6
            
            # Test canonical 4-rotation precession angles 
            ϵ₀, ψₐ, ωₐ, χₐ = precession_angles_rot4(m, tt_c)

            @test r2a*abs(ϵ₀ - ϵe) ≤ 1e-6
            @test r2a*abs(ψₐ - ψe) ≤ 1e-6
            @test r2a*abs(ωₐ - ωe) ≤ 1e-6
            @test r2a*abs(χₐ - χe) ≤ 1e-6

            # Precession matrix 
            Pe = bp06(DJ2000, tt_d)[2]
            Pi = orient_precession(m, tt_c)
            @test v2as(Pe*v, Pi*v) ≤ 1e-6
                        
        end

    end

    @testset "Fukushima-Williams" verbose=true begin 

        for j in 1:50 

            tt_c = rand()
            j == 1 && (tt_c = 1)
            tt_d = tt_c*Tempo.CENTURY2DAY

            v = rand(BigFloat, 3)
            v /= norm(v)

            # -- Fukushima-Williams angles (< 0.1 μas)
            fw = fw_angles(iers2010a, tt_c)
            fe = pfw06(DJ2000, tt_d)

            for j in 1:4
                @test r2a*abs(fw[j] - fe[j]) ≤ 1e-7
            end

            # -- FW Rotation Matrix (< 0.1 μas)
            Rₑ = fw2m(fe[1], fe[2], fe[3], fe[4])
            Rₐ = fw_matrix(fw[1], fw[2], fw[3], fw[4])

            @test v2as(Rₑ * v, Rₐ * v) ≤ 1e-7

            # -- FW to (X, Y) CIP coordinates (< 0.1 μas)
            x, y = fw2xy(fw[1], fw[2], fw[3], fw[4])
            xe, ye = ERFA.fw2xy(fe[1], fe[2], fe[3], fe[4])

            @test r2a*abs(xe-x) ≤ 1e-7
            @test r2a*abs(ye-y) ≤ 1e-7

        end

    end

    @testset "Bias-Precession-Nutation" verbose=true begin 

        for _ in 1:50

            tt_c = rand()/4
            tt_d = tt_c*Tempo.CENTURY2DAY

            v = rand(BigFloat, 3)
            v /= norm(v)

            # --- Testing IERS 1996 models

            # Test frame bias matrix (should be the identity)
            @test orient_bias(iers1996, tt_c) == DCM(1I)

            # Test bias-precession matrix (which equals the precession matrix)
            Pe = pmat76(DJ2000, tt_d)
            PB = orient_pb(iers1996, tt_c)
            @test v2as(Pe*v, PB*v) ≤ 1e-6

            # Test nutation-precession-bias matrix 
            PNe = pnm80(DJ2000, tt_d)
            NPB = orient_npb(iers1996, tt_c)
            @test v2as(PNe*v, NPB*v) ≤ 1e-6

            # --- Testing IERS 2003A/B models 
            Be, _, PBe = bp00(DJ2000, tt_d)  

            # Test common functions
            for m in (iers2003a, iers2003b)
                # Bias matrix (< 1 μas)
                B = orient_bias(m, tt_c) 
                @test v2as(Be*v, B*v) ≤ 1e-6            

                # Precession-Bias matrix (< 1 μas)
                PB = orient_pb(m, tt_c)
                @test v2as(PBe*v, PB*v) ≤ 1e-6
            end

            # Test Nutation-Precession-Bias matrix 
            NPBe = pnm00a(DJ2000, tt_d)

            # 2003A model < 1 μas
            NPB = orient_npb(iers2003a, tt_c)
            @test v2as(NPBe*v, NPB*v) ≤ 1e-6

            # 2003B model < 2.5 mas
            NPB = orient_npb(iers2003b, tt_c)
            @test v2as(NPBe*v, NPB*v) ≤ 2.5e-3

            # --- Testing IERS 2010 models 
            Be, _, PBe = bp06(DJ2000, tt_d)

            # Test common functions 
            for m in (iers2010a, iers2010b, CPNc, CPNd)
                # Bias matrix (< 1 μas)
                B = orient_bias(m, tt_c) 
                @test v2as(Be*v, B*v) ≤ 1e-6            

                # Precession-Bias matrix (< 1 μas)
                PB = orient_pb(m, tt_c)
                @test v2as(PBe*v, PB*v) ≤ 1e-6
            end

            # Test Nutation-Precession-Bias matrix 
            NPBe = pnm06a(DJ2000, tt_d)

            # 2010A model < 1 μas
            NPB = orient_npb(iers2010a, tt_c)
            @test v2as(NPBe*v, NPB*v) ≤ 1e-6

            # 2010B model < 2.5 mas
            NPB = orient_npb(iers2010b, tt_c)
            @test v2as(NPBe*v, NPB*v) ≤ 2.5e-3

            # CPNc model < 30 mas
            NPB = orient_npb(CPNc, tt_c)
            @test v2as(NPBe*v, NPB*v) ≤ 30e-3

            # CPNd model < 1 mas
            NPB = orient_npb(CPNd, tt_c)
            @test v2as(NPBe*v, NPB*v) ≤ 1
        end

    end

    @testset "CIP" verbose=true begin 

        for _ in 1:50

            tt_c = rand()/2
            tt_d = tt_c*Tempo.CENTURY2DAY

            v = rand(BigFloat, 3)
            v /= norm(v)

            # --- Testing IERS 1996 (< 0.1 mas)
            x, y = cip_xy(iers1996, tt_c)
            xe, ye = ERFA.bpn2xy(pnm80(DJ2000, tt_d))
            
            @test r2a*abs(xe-x) ≤ 1e-4 
            @test r2a*abs(ye-y) ≤ 1e-4
            
            # TODO: test on cio locator is missing

            # CIP vector 
            cip = cip_vector(iers1996, tt_c)
            @test r2a*abs(cip[1] - xe)      ≤ 15e-3
            @test r2a*abs(cip[2] - ye)      ≤ 15e-3
            @test r2a*abs(cip[3] - sqrt(1 - xe^2 + ye^2)) ≤ 15e-3

            # --- Testing IERS 2003A model (< 1 μas)
            # CIP coordinates
            x, y = cip_xy(iers2003a, tt_c)
            xe, ye = bpn2xy(pnm00a(DJ2000, tt_d))

            @test r2a*abs(xe-x) ≤ 1e-7
            @test r2a*abs(ye-y) ≤ 1e-7

            # CIO coordinates 
            s = cip_xys(iers2003a, tt_c)[3]
            se = s00a(DJ2000, tt_d)
            @test r2a*abs(se-s) ≤ 1e-7

            # CIP motion matrix
            Qe = c2i00a(DJ2000, tt_d)
            Q = orient_cip_motion(iers2003a, tt_c)
            @test v2as(Qe*v, Q*v) ≤ 1e-7

            # CIP vector 
            cip = cip_vector(iers2003a, tt_c)
            @test r2a*abs(cip[1] - xe)      ≤ 1e-7
            @test r2a*abs(cip[2] - ye)      ≤ 1e-7
            @test r2a*abs(cip[3] - Q[3, 3]) ≤ 1e-7


            # --- Testing IERS 2003B model (< 2.5 mas)
            # CIP coordinates
            x, y = cip_xy(iers2003b, tt_c)

            @test r2a*abs(xe-x) ≤ 2.5e-3
            @test r2a*abs(ye-y) ≤ 2.5e-3

            # CIO coordinates 
            s = cip_xys(iers2003b, tt_c)[3]
            @test r2a*abs(se-s) ≤ 2.5e-3

            # CIP motion matrix
            Qe = c2i00a(DJ2000, tt_d)
            @test v2as(Qe*v, Q*v) ≤ 2.5e-3

            # CIP vector 
            cip = cip_vector(iers2003b, tt_c)
            @test r2a*abs(cip[1] - xe)      ≤ 2.5e-3
            @test r2a*abs(cip[2] - ye)      ≤ 2.5e-3
            @test r2a*abs(cip[3] - Q[3, 3]) ≤ 2.5e-3


            # --- Testing IERS 2010A model (< 1 μas)
            # CIP coordinates
            x, y = cip_xy(iers2010a, tt_c)
            xe, ye = bpn2xy(pnm06a(DJ2000, tt_d))

            @test r2a*abs(xe-x) ≤ 1e-7
            @test r2a*abs(ye-y) ≤ 1e-7

            # CIO coordinates 
            s = cip_xys(iers2010a, tt_c)[3]
            se = s06a(DJ2000, tt_d)
            @test r2a*abs(se-s) ≤ 1e-7

            # CIP motion matrix
            Qe = c2i06a(DJ2000, tt_d)
            Q = orient_cip_motion(iers2010a, tt_c)
            @test v2as(Qe*v, Q*v) ≤ 1e-7

            # CIP vector 
            cip = cip_vector(iers2010a, tt_c)
            @test r2a*abs(cip[1] - xe)      ≤ 1e-7
            @test r2a*abs(cip[2] - ye)      ≤ 1e-7
            @test r2a*abs(cip[3] - Q[3, 3]) ≤ 1e-7


            # --- Testing IERS 2010B model (< 2.5 mas)
            # CIP coordinates
            x, y = cip_xy(iers2010b, tt_c)

            @test r2a*abs(xe-x) ≤ 2.5e-3
            @test r2a*abs(ye-y) ≤ 2.5e-3

            # CIO coordinates 
            s = cip_xys(iers2010b, tt_c)[3]
            @test r2a*abs(se-s) ≤ 2.5e-3

            # CIP motion matrix
            Q = orient_cip_motion(iers2010b, tt_c)
            @test v2as(Qe*v, Q*v) ≤ 2.5e-3

            # CIP vector 
            cip = cip_vector(iers2010b, tt_c)
            @test r2a*abs(cip[1] - xe)      ≤ 2.5e-3
            @test r2a*abs(cip[2] - ye)      ≤ 2.5e-3
            @test r2a*abs(cip[3] - Q[3, 3]) ≤ 2.5e-3


            # --- Testing CPNc model (< 15 mas)
            # CIP coordinates
            x, y = cip_xy(CPNc, tt_c)

            @test r2a*abs(xe-x) ≤ 15e-3
            @test r2a*abs(ye-y) ≤ 15e-3

            # CIO coordinates 
            s = cip_xys(CPNc, tt_c)[3]
            @test r2a*abs(se-s) ≤ 15e-3

            # CIP motion matrix
            Q = orient_cip_motion(CPNc, tt_c)
            @test v2as(Qe*v, Q*v) ≤ 15e-3

            # CIP vector 
            cip = cip_vector(CPNc, tt_c)
            @test r2a*abs(cip[1] - xe)      ≤ 15e-3
            @test r2a*abs(cip[2] - ye)      ≤ 15e-3
            @test r2a*abs(cip[3] - Q[3, 3]) ≤ 15e-3
            

            # --- Testing CPNd model (< 1 as)
            # CIP coordinates
            x, y = cip_xy(CPNd, tt_c)

            @test r2a*abs(xe-x) ≤ 1
            @test r2a*abs(ye-y) ≤ 1

            # CIO coordinates 
            s = cip_xys(CPNd, tt_c)[3]
            @test r2a*abs(se-s) ≤ 1

            # CIP motion matrix
            Q = orient_cip_motion(CPNd, tt_c)
            @test v2as(Qe*v, Q*v) ≤ 1

            # CIP vector 
            cip = cip_vector(CPNd, tt_c)
            @test r2a*abs(cip[1] - xe)      ≤ 1
            @test r2a*abs(cip[2] - ye)      ≤ 1
            @test r2a*abs(cip[3] - Q[3, 3]) ≤ 1
            
        end

    end

    @testset "ERA" verbose=true begin 

        # Test on ERA angles (< 1 μas)
        for _ in 1:50

            tt_c = rand()/2
            tt_d = tt_c*Tempo.CENTURY2DAY

            ERA = orient_era(iers2010a, tt_d)
            ERAₑ = era00(DJ2000, tt_d)

            @test r2a*abs(ERA - ERAₑ) ≤ 1e-7
        end

    end

    @testset "Sidereal Time" verbose=true begin 

        for _ in 1:50 

            # TODO: for these tests we currently assume that the 
            # offset between UTC and UT1 is null! It will have to be 
            # changed once we thrown in the actual EOP values!

            tt_c = rand()/4
            tt_d = tt_c*Tempo.CENTURY2DAY

            ut1_d = tt_d + offset_tt2ut1(tt_d*Tempo.DAY2SEC)/Tempo.DAY2SEC

            # --- Testing 1996 model (< 10 μas)
            # GMST 
            gm  = orient_gmst(iers1996, tt_c)
            gme = gmst82(DJ2000, ut1_d)
            @test r2a*abs(gm-gme) ≤ 1e-5

            # Equation of the equinoxes 
            eeq = equation_equinoxes(iers1996, tt_c)
            eeqe = eqeq94(DJ2000, tt_d)
            @test r2a*abs(eeq-eeqe) ≤ 1e-7

            # GAST (the lower tolerances in the GAST is due to the fact that 
            # ERFA neglects the differences between UT1 and TDB in the 
            # computation of the equation of the equinoxes. 
            ga  = orient_gast(iers1996, tt_c)
            gst = gst94(DJ2000, ut1_d)
            @test r2a*abs(ga-gst) ≤ 1e-3


            # --- Testing 2003A model (< 1 μas)
            # GMST 
            gm  = orient_gmst(iers2003a, tt_c)
            gme = gmst00(DJ2000, ut1_d, DJ2000, tt_d)
            @test r2a*abs(gm-gme) ≤ 1e-6

            # Equation of the equinoxes 
            eeq = equation_equinoxes(iers2003a, tt_c)
            eeqe = ee00a(DJ2000, tt_d)
            @test r2a*abs(eeq-eeqe) ≤ 1e-7

            # GAST
            ga  = orient_gast(iers2003a, tt_c)
            gst = gst00a(DJ2000, ut1_d, DJ2000, tt_d)
            @test r2a*abs(ga-gst) ≤ 1e-7


            # --- Testing 2003B model (< 2.5 mas)
            # GMST 
            gm  = orient_gmst(iers2003b, tt_c)
            @test r2a*abs(gm-gme) ≤ 1e-6

            # Equation of the equinoxes 
            eeq = equation_equinoxes(iers2003b, tt_c)
            @test r2a*abs(eeq-eeqe) ≤ 2.5e-3

            # GAST
            ga  = orient_gast(iers2003b, tt_c)
            @test r2a*abs(ga-gst) ≤ 2.5e-3


            # --- Testing 2010A model (< 1 μas)
            # GMST 
            gm  = orient_gmst(iers2010a, tt_c)
            gme = gmst06(DJ2000, ut1_d, DJ2000, tt_d)
            @test r2a*abs(gm-gme) ≤ 1e-6

            # GAST
            ga  = orient_gast(iers2010a, tt_c)
            gst = gst06a(DJ2000, ut1_d, DJ2000, tt_d)
            @test r2a*abs(ga-gst) ≤ 1e-6


            # --- Testing 2010B model (< 2.5 mas)
            # GMST 
            gm  = orient_gmst(iers2010b, tt_c)
            @test r2a*abs(gm-gme) ≤ 2.5e-3

            # GAST
            ga  = orient_gast(iers2010b, tt_c)
            @test r2a*abs(ga-gst) ≤ 2.5e-3


            # --- Testing CPNc model (< 25 mas)
            # GMST 
            gm  = orient_gmst(CPNc, tt_c)
            @test r2a*abs(gm-gme) ≤ 25e-3

            # GAST
            ga  = orient_gast(CPNc, tt_c)
            @test r2a*abs(ga-gst) ≤ 25e-3


            # --- Testing CPNd model (< 15 mas)
            # GMST 
            gm  = orient_gmst(CPNd, tt_c)
            @test r2a*abs(gm-gme) ≤ 1

            # GAST
            ga  = orient_gast(CPNd, tt_c)
            @test r2a*abs(ga-gst) ≤ 1

        end

    end

    @testset "Polar Motion" verbose=true begin 

        for _ in 1:50

            m = rand((iers1996, iers2003a, iers2003b, iers2010a, iers2010b, CPNc, CPNd))

            tt_c = rand()
            tt_d = tt_c*Tempo.CENTURY2DAY

            v = rand(BigFloat, 3)

            # -- Testing TIO Locator 
            sp = tio_locator(iers2010a, tt_c) 
            spₑ = sp00(DJ2000, tt_d) 

            @test r2a*abs(sp-spₑ) ≤ 1e-6

            @test tio_locator(iers1996, tt_c) == 0
            @test tio_locator(CPNd, tt_c) == 0

            # -- Testing Polar Motion 
            xₚ, yₚ = rand(), rand()

            RPₐ = orient_polar_motion(iers2010a, xₚ, yₚ, tt_c)
            RPₑ = pom00(xₚ, yₚ, sp00(DJ2000, tt_d))

            @test v2as(RPₑ * v, RPₐ * v) ≤ 1e-8

        end
    end

    @testset "Full Rotation without EOP" verbose=true begin 
        for _ in 1:100 

            tt_c = rand() 
            tt_d = tt_c*Tempo.CENTURY2DAY
            
            v = rand(BigFloat, 3)
            v /= norm(v)

            # TODO: write these tests 

        end
    end

    @testset "Full Rotation with EOP" verbose=true begin 
        # TODO: write these tests 

    end

end;

# Function to compute the angle between 2 vectors in arcseconds 
v2as = (x, y) -> acosd(max(-1, min(1, dot(x / norm(x), y / norm(y))))) * 3600


ep = Epoch("2007-04-05T00:00:00 UTC")
ep_tt = convert(TT, ep)

tt_c = j2000c(ep_tt)
tt_d = j2000(ep_tt)

ut1_d = tt_d + offset_tt2ut1(j2000s(ep_tt))/Tempo.DAY2SEC

m = iers1996

begin 

    # CIO-based rotations 
    Q = orient_cip_motion(m, tt_c)
    R = orient_era_rotm(m, ut1_d)

    r_GCRF = randn(BigFloat, 3)
    r_GCRF /= norm(r_GCRF)

    r_CIRF = Q*r_GCRF
    r_TIRF = R*r_CIRF  

    # Equinox-based rotations 
    NPB = orient_npb(m, tt_c)
    R   = angle_to_dcm(orient_gast(m, tt_c), :Z)

    r_TOD  = NPB*r_GCRF
    r_GTOD = R*r_TOD 

    v2as(r_TIRF, r_GTOD)

end

# ITRF-to-PEF 
r_ITRF = [-1033.4793830, 7901.2952754, 6380.3565958]

xₚ = -arcsec2rad(0.140682)
yₚ =  arcsec2rad(0.333309)

W1 = fcn1(xₚ, yₚ)
W2 = fcn2(xₚ, yₚ)
W3 = fcn3(xₚ, yₚ)

a = W1*r_ITRF
b = W2*r_ITRF
c = W3*r_ITRF

v2as = (x, y) -> acosd(max(-1, min(1, dot(x / norm(x), y / norm(y))))) * 3600

t = 0:0.001:1
n = length(t)

err = zeros(n);

v = randn(BigFloat, 3)

for j = 1:n 

    s = tio_locator(iers2010a, t[j])
    R = angle_to_dcm(s, :Z) 

    err[j] = v2as(v, R*v)
end

using GLMakie