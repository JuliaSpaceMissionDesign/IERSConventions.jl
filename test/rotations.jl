
test_dir = artifact"testdata"

# Common routines
r2a = 180 / π * 3600
v2as = (x, y) -> acosd(max(-1, min(1, dot(x / norm(x), y / norm(y))))) * 3600


@testset "Rotations" verbose=true begin 

    # Load EOP data
    eop_load_data!(iers2010a, joinpath(@__DIR__, "assets", "eopc04_20.1962-now.eop.dat"))

    @testset "CIO-based approach" verbose=true begin 
        @testset "CIRF-to-GCRF" verbose=true begin 

            # Retrieve the data
            data = readdlm(joinpath(test_dir, "obspm-cio", "obspm-cipmotion.txt"); skipstart=2)
            n = length(axes(data, 1))

            for j = 1:n
                row = data[j, :]

                ep_utc = Epoch(
                        "$(Int(row[2]))-$(Int(row[3]))-$(Int(row[4]))T"*
                        "$(Int(row[5])):$(Int(row[6])):$(row[7]) UTC"
                )

                ep_tt = convert(TT, ep_utc)

                # The transpose is needed to get the matrix in row-major! 
                R1 = DCM(row[8:16])'
                R2 = iers_rot3_gcrf_to_cirf(j2000s(ep_tt), iers2010a)'

                for _ in 1:10
                    v = rand(BigFloat, 3)
                    @test v2as(R1*v, R2*v) ≤ 3e-6
                end

            end

        end

        @testset "TIRF-to-CIRF without ΔUT1" verbose=true begin 
            data = readdlm(joinpath(test_dir, "obspm-cio", "obspm-era-nodut1.txt"); skipstart=2)
            n = length(axes(data, 1))

            for j = 1:n 
                row = data[j, :]

                ep_utc = Epoch(
                        "$(Int(row[2]))-$(Int(row[3]))-$(Int(row[4]))T"*
                        "$(Int(row[5])):$(Int(row[6])):$(row[7]) UTC"
                )

                ep_ut1 = convert(UT1, ep_utc)

                R1 = DCM(row[8:16])'
                R2 = iers_era_rotm(iers2010a, j2000(ep_utc))'

                for _ in 1:10
                    v = rand(BigFloat, 3)
                    @test v2as(R1*v, R2*v) ≤ 1e-6 
                end
            end
        end

        @testset "TIRF-to-CIRF with ΔUT1" verbose=true begin 

            # Retrieve the data 
            data = readdlm(joinpath(test_dir, "obspm-cio", "obspm-era.txt"); skipstart=2)
            n = length(axes(data, 1))

            for j = 1:n 
                row = data[j, :]

                ep_utc = Epoch(
                        "$(Int(row[2]))-$(Int(row[3]))-$(Int(row[4]))T"*
                        "$(Int(row[5])):$(Int(row[6])):$(row[7]) UTC"
                )

                ep_ut1 = convert(UT1, ep_utc)

                R1 = DCM(row[8:16])'
                R2 = iers_era_rotm(iers2010a, j2000(ep_ut1))'

                for _ in 1:10
                    v = rand(BigFloat, 3)
                    @test v2as(R1*v, R2*v) ≤ 20e-6 
                end
            end

        end

        @testset "ITRF-to-TIRF" verbose=true begin 

            # Retrieve the data
            data = readdlm(joinpath(test_dir, "obspm-cio", "obspm-pm.txt"); skipstart=2)
            n = length(axes(data, 1))
        
            for j = 1:n
                row = data[j, :]

                ep_utc = Epoch(
                        "$(Int(row[2]))-$(Int(row[3]))-$(Int(row[4]))T"*
                        "$(Int(row[5])):$(Int(row[6])):$(row[7]) UTC"
                )

                ep_tt = convert(TT, ep_utc)

                R1 = DCM(row[8:16])'
                R2 = iers_rot3_itrf_to_tirf(j2000s(ep_tt), iers2010a)

                for _ in 1:10
                    v = rand(BigFloat, 3)
                    @test v2as(R1*v, R2*v) ≤ 3e-6
                end

            end
        
        end

        @testset "ITRF-to-CIRF" verbose=true begin 

            # Retrieve the data
            data = readdlm(joinpath(test_dir, "obspm-cio", "obspm-erapm.txt"); skipstart=2)
            n = length(axes(data, 1))
        
            for j = 1:n
                row = data[j, :]

                ep_utc = Epoch(
                        "$(Int(row[2]))-$(Int(row[3]))-$(Int(row[4]))T"*
                        "$(Int(row[5])):$(Int(row[6])):$(row[7]) UTC"
                )

                ep_tt = convert(TT, ep_utc)

                R1 = DCM(row[8:16])'
                R2 = iers_rot3_itrf_to_cirf(j2000s(ep_tt), iers2010a)

                for _ in 1:10
                    v = rand(BigFloat, 3)
                    @test v2as(R1*v, R2*v) ≤ 20e-6
                end

            end
        
        end

        @testset "GCRF-to-ITRF" verbose=true begin 
            # Retrieve the data
            data = readdlm(joinpath(test_dir, "obspm-cio", "obspm-full.txt"); skipstart=2)
            n = length(axes(data, 1))
        
            for j = 1:n
                row = data[j, :]

                ep_utc = Epoch(
                        "$(Int(row[2]))-$(Int(row[3]))-$(Int(row[4]))T"*
                        "$(Int(row[5])):$(Int(row[6])):$(row[7]) UTC"
                )

                ep_tt = convert(TT, ep_utc)

                R  = DCM(row[8:16])'
                Ra = iers_rot3_gcrf_to_itrf(j2000s(ep_tt), iers2010a)'

                # Once we guarantee that the IAU 2006/2000 model is correct, we test all the 
                # other 2010 approximations against that (within their accuracy)
                Rb = iers_rot3_gcrf_to_itrf(j2000s(ep_tt), iers2010b)'
                Rc = iers_rot3_gcrf_to_itrf(j2000s(ep_tt), CPNc)'
                Rd = iers_rot3_gcrf_to_itrf(j2000s(ep_tt), CPNd)' 

                for _ in 1:10
                    v = rand(BigFloat, 3)
                    @test v2as(R*v, Ra*v) ≤ 20e-6
                    @test v2as(Ra*v, Rb*v) ≤ 3e-3
                    @test v2as(Ra*v, Rc*v) ≤ 50e-3
                    @test v2as(Ra*v, Rd*v) ≤ 2

                end

            end

        end

    end 

    @testset "Equinox-based approach" verbose=true begin 

        ep1 = convert(TT, Epoch("2000-01-01T00:00:00 UTC"))
        ep2 = convert(TT, Epoch("2020-01-01T00:00:00 UTC"))
        epochs = ep1:86400:ep2

        # To test the next rotations, we are using data retrieved from the matrix 
        # calculator available on the USNO website, which allows to retrieve the 
        # rotations for the TN36 Equinox-based approach.

        # The EOP data is retrieved from the finals2000A file available at this epoch:
        # 2024-01-09T21:00:00 UTC 
        eop_load_data!(iers2010a, joinpath(@__DIR__, "assets", "finals2000A.data.eop.dat"))

        # Since the USNO matrix calculator does not return matrices with fixed absolute 
        # accuracy (it always displays the coefficients in scientific notation with 9 
        # decimals), we are not able to test the rotation accuracy up to 1μas tolerances 
        # for all epochs. The epoch here select minimises the value of GAST, thus the small 
        # values of the coefficients allow to test the rotations to higher absolute 
        # tolerances! 

        ep_utc = Epoch("2009-09-21T00:00:00 UTC")
        ep_tt = convert(TT, ep_utc)

        # Bias matrix (from EME2000 to ICRF)
        B = [ 1.00000000e+00  7.07836869e-08 -8.05621421e-08;
             -7.07836896e-08  1.00000000e+00 -3.30594317e-08;
              8.05621398e-08  3.30594374e-08  1.00000000e+00];
            
        # Precession matrix (from MOD to EME2000)
        P = [ 9.99997192e-01  2.17365687e-03  9.44504641e-04;
             -2.17365686e-03  9.99997638e-01 -1.03863700e-06;
             -9.44504667e-04 -1.01439491e-06  9.99999554e-01];

        # Bias-Precession-Nutation (from TOD to ICRF)
        BPN = [ 9.99997017e-01  2.24017957e-03  9.73295316e-04; 
               -2.24020240e-03  9.99997490e-01  2.23662098e-05;
               -9.73242769e-04 -2.45465216e-05  9.99999526e-01];

        # GAST matrix (from GTOD to TOD)
        GA = [9.99999994e-01 -1.09254238e-04  0.00000000e+00;
              1.09254238e-04  9.99999994e-01  0.00000000e+00;
              0.00000000e+00  0.00000000e+00  1.00000000e+00];
    
        @testset "GCRF-to-MOD" verbose=true begin
        
            M = iers_rot3_gcrf_to_mod(j2000s(ep_tt), iers2010a)
            Mex = (B*P)'

            @test abs(Mex[1, 1] - M[1, 1]) ≤ 1e-9
            @test abs(Mex[2, 2] - M[2, 2]) ≤ 1e-9 
            @test abs(Mex[3, 3] - M[3, 3]) ≤ 1e-10

            @test abs(Mex[1, 2] - M[1, 2]) ≤ 1e-12
            @test abs(Mex[1, 3] - M[1, 3]) ≤ 1e-12
            @test abs(Mex[2, 1] - M[2, 1]) ≤ 1e-12
            @test abs(Mex[2, 3] - M[2, 3]) ≤ 1e-12
            @test abs(Mex[3, 2] - M[3, 2]) ≤ 1e-12

        end 

        @testset "GCRF-to-TOD" verbose=true begin 
            
            T = iers_rot3_gcrf_to_tod(j2000s(ep_tt), iers2010a)
            Tex = BPN'

            @test abs(Tex[1, 1] - T[1, 1]) ≤ 1e-9
            @test abs(Tex[2, 2] - T[2, 2]) ≤ 1e-9 
            @test abs(Tex[3, 3] - T[3, 3]) ≤ 1e-10

            @test abs(Tex[1, 2] - T[1, 2]) ≤ 1e-12
            @test abs(Tex[1, 3] - T[1, 3]) ≤ 1e-12
            @test abs(Tex[2, 1] - T[2, 1]) ≤ 1e-12
            @test abs(Tex[2, 3] - T[2, 3]) ≤ 1e-12
            @test abs(Tex[3, 2] - T[3, 2]) ≤ 1e-12
            
            # We now also test that with other model approximations, the errors 
            # are within the expected model accuracy! 

            for j in eachindex(epochs)
        
                ep = epochs[j]
                tt_s = j2000s(ep)
        
                S1 = iers_rot3_gcrf_to_tod(tt_s, iers2010a)
        
                v = rand(BigFloat, 3)

                S2 = iers_rot3_gcrf_to_tod(tt_s, iers2010b)
                S3 = iers_rot3_gcrf_to_tod(tt_s, CPNc)
                S4 = iers_rot3_gcrf_to_tod(tt_s, CPNd)
    
                @test v2as(S1*v, S2*v) ≤ 3e-3
                @test v2as(S1*v, S3*v) ≤ 50e-3
                @test v2as(S1*v, S4*v) ≤ 1

            end


        end

        @testset "GCRF-to-GTOD" verbose=true begin 

            G = iers_rot3_gcrf_to_gtod(j2000s(ep_tt), iers2010a)
            Gex = GA'*BPN'

            @test abs(Gex[1, 1] - G[1, 1]) ≤ 1e-9
            @test abs(Gex[2, 2] - G[2, 2]) ≤ 1e-9 
            @test abs(Gex[3, 3] - G[3, 3]) ≤ 1e-10

            @test abs(Gex[1, 2] - G[1, 2]) ≤ 1e-11
            @test abs(Gex[1, 3] - G[1, 3]) ≤ 1e-12
            @test abs(Gex[2, 1] - G[2, 1]) ≤ 1e-11
            @test abs(Gex[2, 3] - G[2, 3]) ≤ 1e-12
            @test abs(Gex[3, 2] - G[3, 2]) ≤ 1e-12

            # We now test that the GTOD frame obtained equals the TIRF for the IAU 2010A 
            # model. Since we've already validated the GCRF-to-TIRF routine, we don't 
            # need other external test data for this 
            for j in eachindex(epochs)

                tt_s = j2000s(epochs[j])
        
                S1 = iers_rot3_gcrf_to_tirf(tt_s, iers2010a)
                S2 = iers_rot3_gcrf_to_gtod(tt_s, iers2010a)
        
                v = rand(BigFloat, 3)
                @test v2as(S1*v, S2*v) ≤ 10e-6 

                # We now also test that with other model approximations, the errors 
                # are within the expected model accuracy! 
                S3 = iers_rot3_gcrf_to_gtod(tt_s, iers2010b)
                S4 = iers_rot3_gcrf_to_gtod(tt_s, CPNc)
                S5 = iers_rot3_gcrf_to_gtod(tt_s, CPNd)

                @test v2as(S2*v, S3*v) ≤ 3e-3
                @test v2as(S2*v, S4*v) ≤ 50e-3
                @test v2as(S2*v, S5*v) ≤ 1

            end

        end


        # Now that we have evaluated the rotations from the GCRF, we can indirectly 
        # test the ones that start from the ITRF at multiple epochs
        # ===========================================================================

        F = DCM{Float64}[]
        for j in eachindex(epochs)
            push!(F, iers_rot3_gcrf_to_itrf(j2000s(epochs[j]), iers2010a))
        end

        @testset "ITRF-to-MOD" verbose=true begin 
            
            for j in eachindex(epochs)

                v_g = rand(BigFloat, 3)
                v_i = F[j]*v_g; v_i /= norm(v_i)

                tt_s = j2000s(epochs[j])

                S1 = iers_rot3_gcrf_to_mod(tt_s, iers2010a)
                S2 = iers_rot3_itrf_to_mod(tt_s, iers2010a)

                @test v2as(S1*v_g, S2*v_i) ≤ 3e-6

            end

        end

        @testset "ITRF-to-TOD" verbose=true begin 

            for j in eachindex(epochs)

                v_g = rand(BigFloat, 3)
                v_i = F[j]*v_g; v_i /= norm(v_i)

                tt_s = j2000s(epochs[j])

                S1 = iers_rot3_gcrf_to_tod(tt_s, iers2010a)
                S2 = iers_rot3_itrf_to_tod(tt_s, iers2010a)

                @test v2as(S1*v_g, S2*v_i) ≤ 3e-6

                # We now also test that with other model approximations, the errors 
                # are within the expected model accuracy! 

                S3 = iers_rot3_itrf_to_tod(tt_s, iers2010b)
                S4 = iers_rot3_itrf_to_tod(tt_s, CPNc)
                S5 = iers_rot3_itrf_to_tod(tt_s, CPNd)

                @test v2as(S2*v_g, S3*v_g) ≤ 3e-3
                @test v2as(S2*v_g, S4*v_g) ≤ 50e-3
                @test v2as(S2*v_g, S5*v_g) ≤ 1

            end

        end

        @testset "ITRF-to-GTOD" verbose=true begin 

            for j in eachindex(epochs)

                v_g = rand(BigFloat, 3)
                v_i = F[j]*v_g; v_i /= norm(v_i)

                tt_s = j2000s(epochs[j])

                S1 = iers_rot3_gcrf_to_gtod(tt_s, iers2010a)
                S2 = iers_rot3_itrf_to_gtod(tt_s, iers2010a)

                @test v2as(S1*v_g, S2*v_i) ≤ 3e-6

                # We now also test that with other model approximations, the errors 
                # are within the expected model accuracy! 

                S3 = iers_rot3_itrf_to_gtod(tt_s, iers2010b)
                S4 = iers_rot3_itrf_to_gtod(tt_s, CPNc)
                S5 = iers_rot3_itrf_to_gtod(tt_s, CPNd)

                @test v2as(S2*v_g, S3*v_g) ≤ 3e-3
                @test v2as(S2*v_g, S4*v_g) ≤ 50e-3
                @test v2as(S2*v_g, S5*v_g) ≤ 1

            end
        end

        # The issue with the PEF is that we do not have an absolute way to test 
        # this reference frame. Thus we assume that one of the two functions is correct 
        # and test its inverse against it.

        # TODO: is it possible to find a software that provides it?
        
        @testset "ITRF-to-PEF" verbose=true begin 

            for j in eachindex(epochs)

                v_g = rand(BigFloat, 3)
                v_i = F[j]*v_g; v_i /= norm(v_i)

                tt_s = j2000s(epochs[j])

                S1 = iers_rot3_gcrf_to_pef(tt_s, iers2010a)
                S2 = iers_rot3_itrf_to_pef(tt_s, iers2010a)

                @test v2as(S1*v_g, S2*v_i) ≤ 3e-6

                # We are also testing that the error in the ITRF-to-PEF with the CPNd 
                # model is within the expected model accuracy (i.e., we test that the 
                # actual magnitude of the rotation is below 1 arcsec)

                S3 = iers_rot3_itrf_to_pef(tt_s, iers2010b)
                S4 = iers_rot3_itrf_to_pef(tt_s, CPNc)
                S5 = iers_rot3_itrf_to_pef(tt_s, CPNd)

                @test v2as(S2*v_i, S3*v_i) ≤ 3e-3
                @test v2as(S2*v_i, S4*v_i) ≤ 50e-3
                @test v2as(S2*v_i, S5*v_i) ≤ 1
            
            end

        end

    end

end;
