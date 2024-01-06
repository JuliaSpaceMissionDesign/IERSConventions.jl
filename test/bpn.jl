
r2a = 180 / π * 3600
v2as = (x, y) -> acosd(max(-1, min(1, dot(x / norm(x), y / norm(y))))) * 3600

@testset "Obliquity" verbose = true begin

    for _ in 1:50

        tt_c = rand()/2
        tt_d = tt_c*Tempo.CENTURY2DAY

        # For the 2003 models we also add the precession-rate corrections 
        ϵ96 = obl80(DJ2000, tt_d)
        ϵ00 = ϵ96 + pr00(DJ2000, tt_d)[2]
        ϵ06 = obl06(DJ2000, tt_d)
        
        # Accurate up to 1nas
        @test r2a*(abs(iers_obliquity(iers1996, tt_c) - ϵ96))     ≤ 1e-9
        @test r2a*(abs(iers_obliquity(iers2003a, tt_c) - ϵ00))    ≤ 1e-9
        @test r2a*(abs(iers_obliquity(iers2010a, tt_c) - ϵ06))    ≤ 1e-9

    end
end;

@testset "Nutation" verbose=true begin 

    @testset "Nutation Components" verbose = true begin
        for _ in 1:50
            tt_c = rand() 
            tt_d = tt_c*Tempo.CENTURY2DAY

            # -- Testing IERS 1996 model (< 0.1μas)
            Δψ, Δϵ = iers_nutation_comp(iers1996, tt_c)
            Δp, Δe = nut80(DJ2000, tt_d)

            @test r2a*abs(Δψ - Δp) ≤ 1e-7
            @test r2a*abs(Δϵ - Δe) ≤ 1e-7

            # -- Testing IERS 2000A model (< 0.1μas)
            Δψ, Δϵ = iers_nutation_comp(iers2003a, tt_c)
            Δp, Δe = nut00a(DJ2000, tt_d)

            @test r2a*abs(Δψ - Δp) ≤ 1e-7
            @test r2a*abs(Δϵ - Δe) ≤ 1e-7

            # -- Testing IERS 2000B model (< 1mas)
            Δψ, Δϵ = iers_nutation_comp(iers2003b, tt_c)
            Δp, Δe = nut00b(DJ2000, tt_d)

            @test r2a*abs(Δψ - Δp) ≤ 1e-7
            @test r2a*abs(Δϵ - Δe) ≤ 1e-7

            # -- Testing IERS 2010A model (< 0.1μas)
            Δψ, Δϵ = iers_nutation_comp(iers2010a, tt_c)
            Δp, Δe = nut06a(DJ2000, tt_d)

            @test r2a*abs(Δψ - Δp) ≤ 1e-7
            @test r2a*abs(Δϵ - Δe) ≤ 1e-7

            # -- Testing IERS 2010B model (< 3mas)
            Δψ, Δϵ = iers_nutation_comp(iers2010b, tt_c)
            @test r2a*abs(Δψ - Δp) ≤ 3e-3
            @test r2a*abs(Δϵ - Δe) ≤ 3e-3
            
            # -- Testing IERS CPNc model (< 50mas)
            Δψ, Δϵ = iers_nutation_comp(CPNc, tt_c)
            @test r2a*abs(Δψ - Δp) ≤ 50e-3
            @test r2a*abs(Δϵ - Δe) ≤ 50e-3

            # -- Testing IERS CPNd model (< 2as)
            Δψ, Δϵ = iers_nutation_comp(CPNd, tt_c)
            @test r2a*abs(Δψ - Δp) ≤ 2
            @test r2a*abs(Δϵ - Δe) ≤ 2
        end
    end

    @testset "Nutation Matrix" verbose=true begin 
        for _ in 1:50
            tt_c = rand() 
            tt_d = tt_c*Tempo.CENTURY2DAY

            v = rand(BigFloat, 3)

            # --- Testing IERS 1996 model (< 0.1μas)
            N1 = iers_nutation(iers1996, tt_c)
            Ne = nutm80(DJ2000, tt_d)
            @test v2as(N1*v, Ne*v) ≤ 1e-7

            # --- Testing IERS 2003 models

            Ne = num00a(DJ2000, tt_d)
            ve = Ne*v; 

            # IERS 2000A model (< 0.1μas)
            N1 = iers_nutation(iers2003a, tt_c)
            @test v2as(N1*v, ve) ≤ 1e-7

            # IERS 2000B model (< 3mas)
            N1 = iers_nutation(iers2003b, tt_c)
            @test v2as(N1*v, ve) ≤ 3e-3

            # --- Testing IERS 2010 models

            Ne = num06a(DJ2000, tt_d)
            ve = Ne*v; 

            # IERS 2010A model (< 0.1μas)
            N1 = iers_nutation(iers2010a, tt_c)
            @test v2as(N1*v, ve) ≤ 1e-7

            # IERS 2010B model (< 3mas)
            N1 = iers_nutation(iers2010b, tt_c)
            @test v2as(N1*v, ve) ≤ 3e-3

            # CPNc model (< 50mas)
            N1 = iers_nutation(CPNc, tt_c)
            @test v2as(N1*v, ve) ≤ 50e-3

            # CPNd model (< 2as)
            N1 = iers_nutation(CPNd, tt_c)
            @test v2as(N1*v, ve) ≤ 2

        end

    end

end;


@testset "Precession" verbose=true begin 

    for j in 1:50

        tt_c = rand()
        
        j == 1 && (tt_c = 1)
        tt_d = tt_c*Tempo.CENTURY2DAY

        v = rand(BigFloat, 3)
        v /= norm(v)

        # -- Testing 1996 models (< 1 μas)
        Pe = pmat76(DJ2000, tt_d)
        v1 = Pe*v; v1 /= norm(v1); 

        # Test traditional parameterization
        zₐ, θₐ, ζₐ = IERSConventions.precession_angles_rot3(iers1996, tt_c)
        ζe, ze, θe = prec76(DJ2000, 0, DJ2000, tt_d)
        
        @test r2a*abs(zₐ - ze) ≤ 1e-6
        @test r2a*abs(θₐ - θe) ≤ 1e-6
        @test r2a*abs(ζₐ - ζe) ≤ 1e-6

        # Test Capitaine parameterization (since we do not have ERFA expressions 
        # that provide this angles, we check that the resulting precession matrix is 
        # the same as that obtained with the original parameterization)
        ϵ₀, ψₐ, ωₐ, χₐ = IERSConventions.precession_angles_rot4(iers1996, tt_c)
        P = angle_to_dcm(χₐ, :Z)*angle_to_dcm(ϵ₀, -ψₐ, -ωₐ, :XZX)

        # This test is accurate up to 0.1 mas, which is anyway below the precision of 
        # the 1977 precession model. 
        @test v2as(v1, P*v) ≤ 1e-4

        # Precession matrix 
        Pi = iers_precession(iers1996, tt_c)
        @test v2as(v1, Pi*v) ≤ 1e-6

        # --- Testing IERS 2003 models (< 1 μas)
        m = rand((iers2003a, iers2003b))

        Pe = bp00(DJ2000, tt_d)[2]
        v1 = Pe*v; v1 /= norm(v1); 

        # Test traditional parameterization
        zₐ, θₐ, ζₐ = IERSConventions.precession_angles_rot3(iers2003a, tt_c)
        P = angle_to_dcm(-ζₐ, θₐ, -zₐ, :ZYZ)
        @test v2as(v1, P*v) ≤ 1e-6

        # Precession matrix 
        # This one also tests the Capitaine parameterization. 
        Pi = iers_precession(m, tt_c)
        @test v2as(v1, Pi*v) ≤ 1e-6

        # --- Testing IERS 2010 models (< 1 μas)
        m = rand((iers2010a, iers2010b, CPNc, CPNd))

        # Retrieve all the angles 
        ϵe, ψe, ωe, χe, ze, ζe, θe = p06e(DJ2000, tt_d)[[1, 2, 3, 9, 10, 11, 12]]

        # Test equatorial precession angles 
        zₐ, θₐ, ζₐ = IERSConventions.precession_angles_rot3(m, tt_c)

        @test r2a*abs(zₐ - ze) ≤ 1e-6
        @test r2a*abs(θₐ - θe) ≤ 1e-6
        @test r2a*abs(ζₐ - ζe) ≤ 1e-6
        
        # Test canonical 4-rotation precession angles 
        ϵ₀, ψₐ, ωₐ, χₐ = IERSConventions.precession_angles_rot4(m, tt_c)

        @test r2a*abs(ϵ₀ - ϵe) ≤ 1e-6
        @test r2a*abs(ψₐ - ψe) ≤ 1e-6
        @test r2a*abs(ωₐ - ωe) ≤ 1e-6
        @test r2a*abs(χₐ - χe) ≤ 1e-6

        # Precession matrix 
        Pe = bp06(DJ2000, tt_d)[2]
        Pi = iers_precession(m, tt_c)
        @test v2as(Pe*v, Pi*v) ≤ 1e-6
                    
    end

end;

@testset "Fukushima-Williams" verbose=true begin 

    for j in 1:50 

        tt_c = rand()
        j == 1 && (tt_c = 1)
        tt_d = tt_c*Tempo.CENTURY2DAY

        v = rand(BigFloat, 3)
        v /= norm(v)

        # -- Fukushima-Williams angles (< 0.1 μas)
        fw = IERSConventions.fw_angles(iers2010a, tt_c)
        fe = pfw06(DJ2000, tt_d)

        for j in 1:4
            @test r2a*abs(fw[j] - fe[j]) ≤ 1e-7
        end

        # -- FW Rotation Matrix (< 0.1 μas)
        Rₑ = fw2m(fe[1], fe[2], fe[3], fe[4])
        Rₐ = IERSConventions.fw_matrix(fw[1], fw[2], fw[3], fw[4])

        @test v2as(Rₑ * v, Rₐ * v) ≤ 1e-7

        # -- FW to (X, Y) CIP coordinates (< 0.1 μas)
        x, y = IERSConventions.fw2xy(fw[1], fw[2], fw[3], fw[4])
        xe, ye = ERFA.fw2xy(fe[1], fe[2], fe[3], fe[4])

        @test r2a*abs(xe-x) ≤ 1e-7
        @test r2a*abs(ye-y) ≤ 1e-7

    end

end;

@testset "Bias-Precession-Nutation" verbose=true begin 

    for _ in 1:50

        tt_c = rand()
        tt_d = tt_c*Tempo.CENTURY2DAY

        v = rand(BigFloat, 3)
        v /= norm(v)

        # --- Testing IERS 1996 models

        # Test frame bias matrix (should be the identity)
        @test iers_bias(iers1996, tt_c) == DCM(1I)

        # Test bias-precession matrix (which equals the precession matrix)
        Pe = pmat76(DJ2000, tt_d)
        PB = iers_pb(iers1996, tt_c)
        @test v2as(Pe*v, PB*v) ≤ 1e-6

        # Test nutation-precession-bias matrix 
        PNe = pnm80(DJ2000, tt_d)
        NPB = iers_npb(iers1996, tt_c)
        @test v2as(PNe*v, NPB*v) ≤ 1e-6

        # --- Testing IERS 2003A/B models 
        Be, _, PBe = bp00(DJ2000, tt_d)  

        # Test common functions
        for m in (iers2003a, iers2003b)
            # Bias matrix (< 1 μas)
            B = iers_bias(m, tt_c) 
            @test v2as(Be*v, B*v) ≤ 1e-6            

            # Precession-Bias matrix (< 1 μas)
            PB = iers_pb(m, tt_c)
            @test v2as(PBe*v, PB*v) ≤ 1e-6
        end

        # Test Nutation-Precession-Bias matrix 
        NPBe = pnm00a(DJ2000, tt_d)

        # 2003A model < 1 μas
        NPB = iers_npb(iers2003a, tt_c)
        @test v2as(NPBe*v, NPB*v) ≤ 1e-6

        # 2003B model < 2.5 mas
        NPB = iers_npb(iers2003b, tt_c)
        @test v2as(NPBe*v, NPB*v) ≤ 2.5e-3

        # --- Testing IERS 2010 models 
        Be, _, PBe = bp06(DJ2000, tt_d)

        # Test common functions 
        for m in (iers2010a, iers2010b, CPNc, CPNd)
            # Bias matrix (< 1 μas)
            B = iers_bias(m, tt_c) 
            @test v2as(Be*v, B*v) ≤ 1e-6            

            # Precession-Bias matrix (< 1 μas)
            PB = iers_pb(m, tt_c)
            @test v2as(PBe*v, PB*v) ≤ 1e-6
        end

        # Test Nutation-Precession-Bias matrix 
        NPBe = pnm06a(DJ2000, tt_d)

        # 2010A model < 1 μas
        NPB = iers_npb(iers2010a, tt_c)
        @test v2as(NPBe*v, NPB*v) ≤ 1e-6

        # 2010B model < 3 mas
        NPB = iers_npb(iers2010b, tt_c)
        @test v2as(NPBe*v, NPB*v) ≤ 3e-3

        # CPNc model < 50 mas
        NPB = iers_npb(CPNc, tt_c)
        @test v2as(NPBe*v, NPB*v) ≤ 50e-3

        # CPNd model < 2 mas
        NPB = iers_npb(CPNd, tt_c)
        @test v2as(NPBe*v, NPB*v) ≤ 2
    end

end;