
r2a = 180 / π * 3600
v2as = (x, y) -> acosd(max(-1, min(1, dot(x / norm(x), y / norm(y))))) * 3600

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
        Q = iers_cip_motion(iers2003a, tt_c)
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
        Q = iers_cip_motion(iers2010a, tt_c)
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
        Q = iers_cip_motion(iers2010b, tt_c)
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
        Q = iers_cip_motion(CPNc, tt_c)
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
        Q = iers_cip_motion(CPNd, tt_c)
        @test v2as(Qe*v, Q*v) ≤ 1

        # CIP vector 
        cip = cip_vector(CPNd, tt_c)
        @test r2a*abs(cip[1] - xe)      ≤ 1
        @test r2a*abs(cip[2] - ye)      ≤ 1
        @test r2a*abs(cip[3] - Q[3, 3]) ≤ 1
        
    end

end;
