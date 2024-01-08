
test_dir = artifact"testdata"

# Common routines
r2a = 180 / π * 3600
v2as = (x, y) -> acosd(max(-1, min(1, dot(x / norm(x), y / norm(y))))) * 3600

# With that epoch lists, retrieve all the required EOP data 
eopfile = joinpath(test_dir, "eopc04_20.1962-now.txt")
eop_generate_from_txt(iers2010a, eopfile, joinpath(@__DIR__, "assets", "eopc04"))
eop_load_data!(joinpath(@__DIR__, "assets", "eopc04.eop.dat"), iers2010a)

@testset "Rotations" verbose=true begin 

    @testset "CIO-based approach" verbose=true begin 
        @testset "CIRF-to-GCRF" verbose=true begin 

            # Retrieve the data
            data = readdlm(joinpath(test_dir, "obspm-cipmotion.txt"); skipstart=2)
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
            data = readdlm(joinpath(test_dir, "obspm-era-nodut1.txt"); skipstart=2)
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
            data = readdlm(joinpath(test_dir, "obspm-era.txt"); skipstart=2)
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
            data = readdlm(joinpath(test_dir, "obspm-pm.txt"); skipstart=2)
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
            data = readdlm(joinpath(test_dir, "obspm-erapm.txt"); skipstart=2)
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
            data = readdlm(joinpath(test_dir, "obspm-full.txt"); skipstart=2)
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
        eps = ep1:86400:ep2
    
        # We now test that the GTOD frame obtained equals the TIRF for the IAU 2010A model.  
        # Since we've already validated the GCRF-to-TIRF routine, we don't require external 
        # test data for this 
        
        @testset "GCRF-to-GTOD" verbose=true begin 
            for j in eachindex(eps)
        
                ep = eps[j]
        
                S1 = iers_rot3_gcrf_to_tirf(j2000s(ep), iers2010a)
                S2 = iers_rot3_gcrf_to_gtod(j2000s(ep), iers2010a)
        
                for _ in 1:10
                    v = rand(BigFloat, 3)
                    @test v2as(S1*v, S2*v) ≤ 10e-6 
                end
            end
        end

    end

end;
