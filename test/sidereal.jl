
r2a = 180 / π * 3600

@testset "Sidereal Time" verbose=true begin 

    for _ in 1:50 

        # TODO: for these tests we currently assume that the 
        # offset between UTC and UT1 is null! It will have to be 
        # changed once we thrown in the actual EOP values!

        tt_c = rand()/4
        tt_d = tt_c*Tempo.CENTURY2DAY

        ut1_d = tt_d + IERS.offset_tt2ut1(tt_d*Tempo.DAY2SEC)/Tempo.DAY2SEC

        # --- Testing 1996 model (< 10 μas)
        # GMST 
        gm  = iers_gmst(iers1996, tt_c)
        gme = gmst82(DJ2000, ut1_d)
        @test r2a*abs(gm-gme) ≤ 1e-5

        # Equation of the equinoxes 
        eeq = IERS.equation_equinoxes(iers1996, tt_c)
        eeqe = eqeq94(DJ2000, tt_d)
        @test r2a*abs(eeq-eeqe) ≤ 1e-7

        # GAST (the lower tolerances in the GAST is due to the fact that 
        # ERFA neglects the differences between UT1 and TDB in the 
        # computation of the equation of the equinoxes. 
        ga  = iers_gast(iers1996, tt_c)
        gst = gst94(DJ2000, ut1_d)
        @test r2a*abs(ga-gst) ≤ 1e-3


        # --- Testing 2003A model (< 1 μas)
        # GMST 
        gm  = iers_gmst(iers2003a, tt_c)
        gme = gmst00(DJ2000, ut1_d, DJ2000, tt_d)
        @test r2a*abs(gm-gme) ≤ 1e-6

        # Equation of the equinoxes 
        eeq = IERS.equation_equinoxes(iers2003a, tt_c)
        eeqe = ee00a(DJ2000, tt_d)
        @test r2a*abs(eeq-eeqe) ≤ 1e-7

        # GAST
        ga  = iers_gast(iers2003a, tt_c)
        gst = gst00a(DJ2000, ut1_d, DJ2000, tt_d)
        @test r2a*abs(ga-gst) ≤ 1e-7


        # --- Testing 2003B model (< 2.5 mas)
        # GMST 
        gm  = iers_gmst(iers2003b, tt_c)
        @test r2a*abs(gm-gme) ≤ 1e-6

        # Equation of the equinoxes 
        eeq = IERS.equation_equinoxes(iers2003b, tt_c)
        @test r2a*abs(eeq-eeqe) ≤ 2.5e-3

        # GAST
        ga  = iers_gast(iers2003b, tt_c)
        @test r2a*abs(ga-gst) ≤ 2.5e-3


        # --- Testing 2010A model (< 1 μas)
        # GMST 
        gm  = iers_gmst(iers2010a, tt_c)
        gme = gmst06(DJ2000, ut1_d, DJ2000, tt_d)
        @test r2a*abs(gm-gme) ≤ 1e-6

        # GAST
        ga  = iers_gast(iers2010a, tt_c)
        gst = gst06a(DJ2000, ut1_d, DJ2000, tt_d)
        @test r2a*abs(ga-gst) ≤ 1e-6


        # --- Testing 2010B model (< 2.5 mas)
        # GMST 
        gm  = iers_gmst(iers2010b, tt_c)
        @test r2a*abs(gm-gme) ≤ 2.5e-3

        # GAST
        ga  = iers_gast(iers2010b, tt_c)
        @test r2a*abs(ga-gst) ≤ 2.5e-3


        # --- Testing CPNc model (< 30 mas)
        # GMST 
        gm  = iers_gmst(CPNc, tt_c)
        @test r2a*abs(gm-gme) ≤ 30e-3

        # GAST
        ga  = iers_gast(CPNc, tt_c)
        @test r2a*abs(ga-gst) ≤ 30e-3


        # --- Testing CPNd model (< 15 mas)
        # GMST 
        gm  = iers_gmst(CPNd, tt_c)
        @test r2a*abs(gm-gme) ≤ 1

        # GAST
        ga  = iers_gast(CPNd, tt_c)
        @test r2a*abs(ga-gst) ≤ 1

    end

end;