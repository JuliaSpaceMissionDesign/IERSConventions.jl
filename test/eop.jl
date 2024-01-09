
test_dir = artifact"testdata"

get_row(data, mjd) = findfirst(x -> x >= mjd, data[:, 1])

@testset "EOP Routines" verbose=true begin 
    @testset "EOP Parsers" verbose=true begin 

        @testset "CSV Files" verbose=true begin 
        
            # Testing of EOP C04 for ITRF 20
            # ====================================================

            fname = "eopc04_20.1962-now"

            src_path = joinpath(test_dir, "eop", "csv", fname*".csv")
            out_path = joinpath(@__DIR__, "assets", fname)

            eop_generate_from_csv(iers2010a, src_path, out_path)
            data1 = readdlm(out_path*".eop.dat")

            # Manually check that the values of a specific MJD match the ones 
            # we expect!

            ep_utc = Epoch("2007-01-01T00:00:00 UTC")
            row = get_row(data1, j2000(ep_utc))

            # check poles
            @test abs(data1[row, 2] - -0.049506) ≤ 1e-6
            @test abs(data1[row, 3] -  0.347383) ≤ 1e-6
        
            # check dut1 and lod 
            @test abs(data1[row, 4] - 0.0376053) ≤ 1e-6
            @test abs(data1[row, 5] - 0.0007636) ≤ 1e-6

            # check CIP corrections 
            @test abs(data1[row, 6] -  0.000145) ≤ 1e-6
            @test abs(data1[row, 7] - -0.000074) ≤ 1e-6

            # Testing of EOP C04 for ITRF 14
            # ====================================================

            fname = "eopc04_14_IAU2000.62-now"

            src_path = joinpath(test_dir, "eop", "csv", fname*".csv")
            out_path = joinpath(@__DIR__, "assets", fname)

            eop_generate_from_csv(iers2010a, src_path, out_path)
            data2 = readdlm(out_path*".eop.dat")

            # Manually check that the values of a specific MJD match the ones 
            # we expect!

            ep_utc = Epoch("2007-01-01T00:00:00 UTC")
            row = get_row(data2, j2000(ep_utc))

            # check poles
            @test abs(data2[row, 2] - -0.049367) ≤ 1e-6
            @test abs(data2[row, 3] -  0.347325) ≤ 1e-6
        
            # check dut1 and lod 
            @test abs(data2[row, 4] - 0.0376484) ≤ 1e-6
            @test abs(data2[row, 5] - 0.0007459) ≤ 1e-6

            # check CIP corrections 
            @test abs(data2[row, 6] -  0.000077) ≤ 1e-6
            @test abs(data2[row, 7] - -0.000057) ≤ 1e-6

            # Testing of EOP C04 with 1980 corrections! 
            # ====================================================
            
            fname = "eopc04_14.62-now"

            src_path = joinpath(test_dir, "eop", "csv", fname*".csv")
            out_path = joinpath(@__DIR__, "assets", fname)

            eop_generate_from_csv(iers2010a, src_path, out_path)
            data3 = readdlm(out_path*".eop.dat")

            # Manually check that the values of a specific MJD match the ones 
            # we expect!

            ep_utc = Epoch("2007-01-01T00:00:00 UTC")
            row = get_row(data3, j2000(ep_utc))

            # check poles
            @test abs(data3[row, 2] - data2[row, 2]) ≤ 1e-6
            @test abs(data3[row, 3] - data2[row, 3]) ≤ 1e-6
        
            # check dut1 and lod 
            @test abs(data3[row, 4] - data2[row, 4]) ≤ 1e-6
            @test abs(data3[row, 5] - data2[row, 5]) ≤ 1e-6

            # check CIP corrections 
            @test abs(data3[row, 8] - -0.058793) ≤ 1e-6
            @test abs(data3[row, 9] - -0.002040) ≤ 1e-6


            # Testing of FINALS IAU 2000
            # ====================================================

            fname = "finals2000A.data"
            src_path = joinpath(test_dir, "eop", "csv", fname*".csv")
            out_path = joinpath(@__DIR__, "assets", fname)

            eop_generate_from_csv(iers2010a, src_path, out_path)
            data4 = readdlm(out_path*".eop.dat")

            # Manually check that the values of a specific MJD match the ones 
            # we expect!

            ep_utc = Epoch("2007-01-01T00:00:00 UTC")
            row = get_row(data4, j2000(ep_utc))

            # check poles
            @test abs(data4[row, 2] - -0.049474) ≤ 1e-6
            @test abs(data4[row, 3] -  0.347382) ≤ 1e-6
        
            # check dut1 and lod 
            @test abs(data4[row, 4] - 0.0375966) ≤ 1e-6
            @test abs(data4[row, 5] - 1e-3*0.7486) ≤ 1e-6

            # check CIP corrections 
            @test abs(data4[row, 6] - 1e-3*0.089) ≤ 1e-6
            @test abs(data4[row, 7] - 1e-3*-0.197) ≤ 1e-6

            
            # Testing of FINALS IAU 1980
            # ====================================================

            fname = "finals.data"
            src_path = joinpath(test_dir, "eop", "csv", fname*".csv")
            out_path = joinpath(@__DIR__, "assets", fname)

            eop_generate_from_csv(iers2010a, src_path, out_path)
            data5 = readdlm(out_path*".eop.dat")

            # Manually check that the values of a specific MJD match the ones 
            # we expect!

            ep_utc = Epoch("2007-01-01T00:00:00 UTC")
            row = get_row(data5, j2000(ep_utc))

            # check poles
            @test abs(data5[row, 2] - data4[row, 2]) ≤ 1e-6
            @test abs(data5[row, 3] - data4[row, 3]) ≤ 1e-6
        
            # check dut1 and lod 
            @test abs(data5[row, 4] - data4[row, 4]) ≤ 1e-6
            @test abs(data5[row, 5] - data4[row, 5]) ≤ 1e-6

            # check nutation corrections 
            @test abs(data5[row, 8] - 1e-3*-58.779) ≤ 1e-6
            @test abs(data5[row, 9] - 1e-3*-2.140) ≤ 1e-6

        end

        @testset "TXT Files" verbose=true begin 
        
            # Testing of EOP C04 for ITRF 20
            # ====================================================

            fname = "eopc04_20.1962-now"

            src_path = joinpath(test_dir, "eop", "txt", fname*".txt")
            out_path = joinpath(@__DIR__, "assets", fname)

            eop_generate_from_txt(iers2010a, src_path, out_path)
            data1 = readdlm(out_path*".eop.dat")

            # Manually check that the values of a specific MJD match the ones 
            # we expect!

            ep_utc = Epoch("2007-01-01T00:00:00 UTC")
            row = get_row(data1, j2000(ep_utc))

            # check poles
            @test abs(data1[row, 2] - -0.049506) ≤ 1e-6
            @test abs(data1[row, 3] -  0.347383) ≤ 1e-6
        
            # check dut1 and lod 
            @test abs(data1[row, 4] - 0.0376053) ≤ 1e-6
            @test abs(data1[row, 5] - 0.0007636) ≤ 1e-6

            # check CIP corrections 
            @test abs(data1[row, 6] -  0.000145) ≤ 1e-6
            @test abs(data1[row, 7] - -0.000074) ≤ 1e-6

            
            # Testing of EOP C04 for ITRF 14
            # ====================================================

            fname = "eopc04_14_IAU2000.62-now"

            src_path = joinpath(test_dir, "eop", "txt", fname*".txt")
            out_path = joinpath(@__DIR__, "assets", fname)

            eop_generate_from_txt(iers2010a, src_path, out_path)
            data2 = readdlm(out_path*".eop.dat")

            # Manually check that the values of a specific MJD match the ones 
            # we expect!

            ep_utc = Epoch("2007-01-01T00:00:00 UTC")
            row = get_row(data2, j2000(ep_utc))

            # check poles
            @test abs(data2[row, 2] - -0.049367) ≤ 1e-6
            @test abs(data2[row, 3] -  0.347325) ≤ 1e-6
        
            # check dut1 and lod 
            @test abs(data2[row, 4] - 0.0376484) ≤ 1e-6
            @test abs(data2[row, 5] - 0.0007459) ≤ 1e-6

            # check CIP corrections 
            @test abs(data2[row, 6] -  0.000077) ≤ 1e-6
            @test abs(data2[row, 7] - -0.000057) ≤ 1e-6


            # Testing of EOP C04 with 1980 corrections! 
            # ====================================================
            
            fname = "eopc04_14.62-now"

            src_path = joinpath(test_dir, "eop", "txt", fname*".txt")
            out_path = joinpath(@__DIR__, "assets", fname)

            eop_generate_from_txt(iers2010a, src_path, out_path)
            data3 = readdlm(out_path*".eop.dat")

            # Manually check that the values of a specific MJD match the ones 
            # we expect!

            ep_utc = Epoch("2007-01-01T00:00:00 UTC")
            row = get_row(data3, j2000(ep_utc))

            # check poles
            @test abs(data3[row, 2] - data2[row, 2]) ≤ 1e-6
            @test abs(data3[row, 3] - data2[row, 3]) ≤ 1e-6
        
            # check dut1 and lod 
            @test abs(data3[row, 4] - data2[row, 4]) ≤ 1e-6
            @test abs(data3[row, 5] - data2[row, 5]) ≤ 1e-6

            # check CIP corrections 
            @test abs(data3[row, 8] - -0.058793) ≤ 1e-6
            @test abs(data3[row, 9] - -0.002040) ≤ 1e-6


            # Testing of FINALS IAU 2000
            # ====================================================

            fname = "finals2000A.data"
            src_path = joinpath(test_dir, "eop", "txt", fname)
            out_path = joinpath(@__DIR__, "assets", fname)

            eop_generate_from_txt(iers2010a, src_path, out_path)
            data4 = readdlm(out_path*".eop.dat")

            # Manually check that the values of a specific MJD match the ones 
            # we expect!

            ep_utc = Epoch("2007-01-01T00:00:00 UTC")
            row = get_row(data4, j2000(ep_utc))

            # check poles
            @test abs(data4[row, 2] - -0.049474) ≤ 1e-6
            @test abs(data4[row, 3] -  0.347382) ≤ 1e-6
        
            # check dut1 and lod 
            @test abs(data4[row, 4] - 0.0375966) ≤ 1e-6
            @test abs(data4[row, 5] - 1e-3*0.7486) ≤ 1e-6

            # check CIP corrections 
            @test abs(data4[row, 6] - 1e-3*0.089) ≤ 1e-6
            @test abs(data4[row, 7] - 1e-3*-0.197) ≤ 1e-6


            # Testing of FINALS IAU 1980
            # ====================================================

            fname = "finals.data"
            src_path = joinpath(test_dir, "eop", "txt", fname)
            out_path = joinpath(@__DIR__, "assets", fname)

            eop_generate_from_txt(iers2010a, src_path, out_path)
            data5 = readdlm(out_path*".eop.dat")

            # Manually check that the values of a specific MJD match the ones 
            # we expect!

            ep_utc = Epoch("2007-01-01T00:00:00 UTC")
            row = get_row(data5, j2000(ep_utc))

            # check poles
            @test abs(data5[row, 2] - data4[row, 2]) ≤ 1e-6
            @test abs(data5[row, 3] - data4[row, 3]) ≤ 1e-6
        
            # check dut1 and lod 
            @test abs(data5[row, 4] - data4[row, 4]) ≤ 1e-6
            @test abs(data5[row, 5] - data4[row, 5]) ≤ 1e-6

            # check nutation corrections 
            @test abs(data5[row, 8] - 1e-3*-58.779) ≤ 1e-6
            @test abs(data5[row, 9] - 1e-3*-2.140) ≤ 1e-6

        end
    end

    @testset "Miscellaneous" verbose=true begin 

        # Test that if no EOP data has been loaded, we get errors! 
        @test_throws ErrorException eop_filename()
        @test_throws ErrorException IERSConventions.eop_δX(iers2010a, 0)

        @test repr(IERSConventions.IERS_EOP_DATA) == "EOPData()\n"
        @test repr(IERSConventions.IERS_EOP) == "EOPInterpolator(init=false)\n"

        dummy_eop = joinpath(test_dir, "eop", "txt", "finals.data")
        @test_throws ArgumentError eop_load_data!(dummy_eop, iers2010a)

        @info "Initialising EOP data"
        eop_file = joinpath(@__DIR__, "assets", "eopc04_20.1962-now.eop.dat")
        eop_load_data!(eop_file, iers2010a)
        
        @test eop_filename() == eop_file 
    end

end;