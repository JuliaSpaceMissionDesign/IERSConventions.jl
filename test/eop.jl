
using IERSConventions: eop_δX, eop_δY, eop_δΔψ, eop_δΔϵ

test_dir = artifact"testdata"

get_row(data, mjd) = findfirst(x -> x >= mjd, data[:, 1])
v2as = (x, y) -> acosd(max(-1, min(1, dot(x / norm(x), y / norm(y))))) * 3600

function gcrf_to_gtod(m, ttc, δΔψ, δΔϵ)

    # Retrieve precession-bias matrix
    PB = iers_pb(m, ttc)

    # Retrieve MOD-to-GTOD rotation 
    RN = IERSConventions.mod_to_gtod3(ttc, m, δΔψ, δΔϵ)
    return RN*PB

end

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

            # Test that when no data is available, the columns are fitted with 0s 
            ep_utc = Epoch("2024-03-15T00:00:00 UTC")
            row = get_row(data4, j2000(ep_utc))
            
            @test all(data4[row, 5:end] .== 0)

            
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

            # Test that when no data is available, the columns are fitted with 0s 
            ep_utc = Epoch("2024-03-15T00:00:00 UTC")
            row = get_row(data5, j2000(ep_utc))
            
            @test all(data5[row, 5:end] .== 0)
            
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

            # Test that when no data is available, the columns are fitted with 0s 
            ep_utc = Epoch("2024-03-15T00:00:00 UTC")
            row = get_row(data4, j2000(ep_utc))
            
            @test all(data4[row, 5:end] .== 0)
            
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

            # Test that when no data is available, the columns are fitted with 0s 
            ep_utc = Epoch("2024-03-15T00:00:00 UTC")
            row = get_row(data5, j2000(ep_utc))
            
            @test all(data5[row, 5:end] .== 0)

        end
    end

    @testset "Miscellaneous" verbose=true begin 

        # Test that if no EOP data has been loaded, we get errors! 
        @test_throws ErrorException eop_filename()
        @test_throws ErrorException IERSConventions.eop_δX(iers2010a, 0)

        @test repr(IERSConventions.IERS_EOP_DATA) == "EOPData()\n"
        @test repr(IERSConventions.IERS_EOP) == "EOPInterpolator(init=false)\n"

        dummy_eop = joinpath(test_dir, "eop", "txt", "finals.data")
        @test_throws ArgumentError eop_load_data!(iers2010a, dummy_eop)

        @info "Initialising EOP data"
        eop_file = joinpath(@__DIR__, "assets", "eopc04_20.1962-now.eop.dat")
        eop_load_data!(iers2010a, eop_file)

        # Test unloading of EOP data 
        eop_unload_data!()

        @test_throws ErrorException eop_filename()
        @test_throws ErrorException IERSConventions.eop_δX(iers2010a, 0)

        eop_load_data!(iers2010a, eop_file)
        @test eop_filename() == eop_file 

        str = "EOPData(filename=\"$eop_file\", from: -10227.5 (UTC) to 8741.5 (UTC))\n"
        @test repr(IERSConventions.IERS_EOP_DATA) == str

        t1 = IERSConventions.IERS_EOP_DATA.cent_TT[1]
        t2 = IERSConventions.IERS_EOP_DATA.cent_TT[end]

        # Check that outside the boundaries all EOP data is zero 
        for m in (iers2010a, iers2010b, CPNc, CPNd, iers2003a, iers2003b, iers1996)
            
            for t in (t1-1e-10, t2+1e-10)

                @test IERSConventions.eop_δΔψ(iers2010a, t) == 0
                @test IERSConventions.eop_δΔϵ(iers2010a, t) == 0

                @test IERSConventions.eop_δX(iers2010a, t) == 0
                @test IERSConventions.eop_δY(iers2010a, t) == 0

                @test IERSConventions.eop_xp(iers2010a, t) == 0
                @test IERSConventions.eop_yp(iers2010a, t) == 0

                @test IERSConventions.offset_tt2ut1(t*Tempo.CENTURY2SEC) == 0

            end

        end

    end

    @testset "EOP Conversion Functions" verbose=true begin 

        # DISCLAIMER: there are few rare cases over the entire 1972-2023 domain in which 
        # some of these tests are slightly violated.

        r2as = 648000/π
        as2r = π/648000

        # We also read the same .CSV from the artifacts to be able to retrieve the 
        # uncertainties over the celestial pole offsets parameters
        csv_dir = joinpath(test_dir, "eop", "csv")
    
        eop10_filename = "eopc04_14_IAU2000.62-now"
        eop96_filename = "eopc04_14.62-now"
        
        # Here we load the EOP data for the 2010 model
        eop10_filepath = joinpath(@__DIR__, "assets", eop10_filename)
        eop96_filepath = joinpath(@__DIR__, "assets", eop96_filename)

        eop_generate_from_csv(iers2010a, joinpath(csv_dir, eop10_filename*".csv"), eop10_filepath)
        eop_generate_from_csv(iers1996,  joinpath(csv_dir, eop96_filename*".csv"), eop96_filepath)
            
        # This is the starting row index (it skips all dates before 2020)
        id_s = 14882
        id_e = 16881

        # Load EOP data matrices
        eop96_data = readdlm(joinpath(csv_dir, "eopc04_14.62-now.csv"), ';'; header=false)[id_s:id_e, :]
        eop10_data = readdlm(joinpath(csv_dir, "eopc04_14_IAU2000.62-now.csv"), ';'; header=false)[id_s:id_e, :]

        # Retrieve 2010 CIP corrections 
        δX_10i = convert(Vector{Float64}, eop10_data[1:end, end-3])*as2r
        δY_10i = convert(Vector{Float64}, eop10_data[1:end, end-1])*as2r

        σ_δX_10i = convert(Vector{Float64}, eop10_data[1:end, end-2])*as2r
        σ_δY_10i = convert(Vector{Float64}, eop10_data[1:end, end])*as2r

        # Retrieve 1980 nutation corrections
        δΔψ_96i = convert(Vector{Float64}, eop96_data[1:end, end-7])*as2r
        δΔϵ_96i = convert(Vector{Float64}, eop96_data[1:end, end-5])*as2r

        σ_δΔψ_96i = convert(Vector{Float64}, eop96_data[1:end, end-6])*as2r
        σ_δΔϵ_96i = convert(Vector{Float64}, eop96_data[1:end, end-4])*as2r

        @testset "From 2010" verbose=true begin

            eop_load_data!(iers2010a, eop10_filepath*".eop.dat")

            for j in eachindex(σ_δΔϵ_96i)

                row = eop10_data[j, :]
                
                # Retrieve the epoch
                ep_utc = Epoch("MJD $(Int(row[1])) UTC")
                ep_tt  = convert(TT, ep_utc)
                ep_ut1 = convert(UT1, ep_tt)
                
                tt_s = j2000s(ep_tt)
                tt_c = j2000c(ep_tt)

                ut1_d = j2000(ep_ut1)

                # Retrieve CIP corrections
                δX_96, δY_96 = eop_δX(iers1996,  tt_c), eop_δY(iers1996, tt_c) 
                δX_03, δY_03 = eop_δX(iers2003a, tt_c), eop_δY(iers2003a, tt_c) 
                δX_10, δY_10 = eop_δX(iers2010a, tt_c), eop_δY(iers2010a, tt_c) 

                # Retrieve nutation corrections 
                δΔψ_96, δΔϵ_96 = eop_δΔψ(iers1996,  tt_c), eop_δΔϵ(iers1996, tt_c)
                δΔψ_03, δΔϵ_03 = eop_δΔψ(iers2003a, tt_c), eop_δΔϵ(iers2003a, tt_c)
                δΔψ_10, δΔϵ_10 = eop_δΔψ(iers2010a, tt_c), eop_δΔϵ(iers2010a, tt_c)

                # Here we test we are interpolating correctly the 2010 values
                @test abs(δX_10 - δX_10i[j])*r2as ≤ 1e-7
                @test abs(δY_10 - δY_10i[j])*r2as ≤ 1e-7

                # For the 2003 and 1996 corrections, we test that the differences between the two 
                # models are below the error caused by the uncertainty of those corrections. 

                # For the 1996 corrections, we obtain residual differences of about 10 μas due to 
                # the neglection of the CIO locator. What we do is that we test that the uncertainty 
                # over the CIP offsets results in errors that are greater than those induced by not 
                # considering the CIO.

                # 3σ - covers 99% of the cases
                dxa = δX_10i[j] + 3σ_δX_10i[j]
                dya = δY_10i[j] + 3σ_δY_10i[j]

                Qr = iers_rot3_gcrf_to_cirf(tt_s, iers2010a)
                Qu = iers_cip_motion(iers2010a, tt_c, dxa, dya)

                Qb = iers_rot3_gcrf_to_cirf(tt_s, iers2003a)
                Qc = iers_rot3_gcrf_to_cirf(tt_s, iers1996)

                # We test that the differences between Q1 and Q2 are smaller than the differences 
                # caused by the uncertainty (over δX_10, δY_10)
                for _ in 1:10
                    v = rand(BigFloat, 3)
                    vr = Qr*v
                    err = v2as(vr, Qu*v)
                    # 1.5 the error because it is just supposed to be a rule-of-thumb test
                    @test v2as(vr, Qb*v) ≤ 1.5err
                    @test v2as(vr, Qc*v) ≤ 1.5err
                end

                # Test the differences between 1996 models (predicted by us vs predicted by IERS). 
                # This means that any residual differences shouldn't matter! 
                @test abs(δΔψ_96 - δΔψ_96i[j]) ≤ 3σ_δΔψ_96i[j]
                @test abs(δΔϵ_96 - δΔϵ_96i[j]) ≤ 3σ_δΔϵ_96i[j]

                # 3σ - covers 99% of the cases
                dpa = δΔψ_96i[j] + 3σ_δΔψ_96i[j]
                dea = δΔϵ_96i[j] + 3σ_δΔϵ_96i[j]
                            
                # First-of-all we test that the differences between the IERS released 1980 corrections 
                # and the ones I am using cause errors that are smaller than those induced by the 
                # uncertainty over those parameters! 
                Ar = gcrf_to_gtod(iers1996, tt_c, δΔψ_96i[j], δΔϵ_96i[j])
                Au = gcrf_to_gtod(iers1996, tt_c, dpa, dea)

                Ab = iers_rot3_gcrf_to_gtod(tt_s, iers1996)
                
                for _ in 1:10 
                    v = rand(BigFloat, 3)
                    vr = Ar*v
                    @test v2as(vr, Ab*v) ≤ v2as(vr, Au*v)
                end

                # Here I could also technically test again that the differences in the 
                # GCRF-to-GTOD matrix caused by the different models are smaller than the differences 
                # induced by the uncertainties over the δX, δY parameters! 
                R = iers_era_rotm(iers2010a, ut1_d)

                Br = R*Qr
                Bu = R*Qu

                Bb = iers_rot3_gcrf_to_gtod(tt_s, iers2003a)
                Bc = iers_rot3_gcrf_to_gtod(tt_s, iers1996)

                for _ in 1:10 
                    v = rand(BigFloat, 3)
                    vr = Br*v
                    err = v2as(vr, Bu*v)

                    @test v2as(vr, Bb*v) ≤ err
                    @test v2as(vr, Bc*v) ≤ err
                end

            end

        end 


    end

end;
