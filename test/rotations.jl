
# Common routines
r2a = 180 / π * 3600
v2as = (x, y) -> acosd(max(-1, min(1, dot(x / norm(x), y / norm(y))))) * 3600

# With that epoch lists, retrieve all the required EOP data 
eopfile = "test/assets/eopc04_20.1972-now.txt"
eop_parse_txt(iers2010a, eopfile, "test/assets/eopc04")
eop_load_data!("test/assets/eopc04.eop.dat", iers2010a)

@testset "Full Rotation with EOP" verbose=true begin 
    # TODO: write these tests 

    @testset "CIP-motion" verbose=true begin 

        # Retrieve the data
        data = readdlm("test/assets/obspm-cipmotion.txt"; skipstart=2)
        n = length(axes(data, 1))

        for j = 1:n
            row = data[j, :]

            ep_utc = Epoch(
                    "$(Int(row[2]))-$(Int(row[3]))-$(Int(row[4]))T"*
                    "$(Int(row[5])):$(Int(row[6])):$(row[7]) UTC"
            )

            ep_tt = convert(TT, ep_utc)

            R1 = DCM(row[8:16])'
            R2 = iers_rot3_gcrf_to_cirf(j2000s(ep_tt), iers2010a)'

            for _ in 1:10
                v = rand(BigFloat, 3)
                @test v2as(R1*v, R2*v) ≤ 3e-6
            end

        end

    end

    @testset "Polar-motion" verbose=true begin 

        # Retrieve the data
        data = readdlm("test/assets/obspm-pm.txt"; skipstart=2)
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

end;
