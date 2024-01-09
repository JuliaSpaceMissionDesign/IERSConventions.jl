
using IERSConventions:  NutationDataParser, 
                        CIODataParser,
                        PoissonSeries,
                        parse_iers_constants, 
                        generate_iers_file

test_dir = artifact"testdata"

@testset "Poisson Series" verbose=true begin 

    # We are testing the data-parser in a sort of unorthodox way. They've already 
    # been used to generate the series coefficients that are available within the 
    # package. Additionally, the functions that are built using those series are 
    # already tested against ERFA functions. So that we know the coefficients constants 
    # are correct. 

    # Therefore, we test the data-parser by verying that if we re-apply them to the 
    # official IERS files, we obtain the same constants we already have within the package.

    @testset "CIO Data Parser" verbose=true begin 
       
        src_file = joinpath(test_dir, "series", "tab5.2d.txt")
        out_file = joinpath(@__DIR__, "assets", "cio2006_tmp.jl")

        isfile(out_file) && rm(out_file)

        parser = CIODataParser(17; fc_del=4, fc_pla=9, c_sin=2, c_cos=3)
    
        data_s = parse_iers_constants(src_file, parser)
        generate_iers_file(out_file, :CIP2006_test, data_s)
        
        include(out_file)
        @test CIP2006_test == IERSConventions.COEFFS_CIO2006_S

    end

    @testset "Nutation Data Parser" verbose=true begin 

        src_file = joinpath(test_dir, "series", "tab5.1.txt")
        out_file = joinpath(@__DIR__, "assets", "nut1980_tmp.jl")

        isfile(out_file) && rm(out_file)

        parser_psi = NutationDataParser(10; fc_del=1, c_sin=[7, 8])
        parser_eps = NutationDataParser(10; fc_del=1, c_cos=[9, 10])

        data_psi = parse_iers_constants(src_file, parser_psi)
        generate_iers_file(out_file, :NUT1980_ψ_test, data_psi)

        data_eps = parse_iers_constants(src_file, parser_eps)
        generate_iers_file(out_file, :NUT1980_ϵ_test, data_eps)

        include(out_file)
        @test NUT1980_ψ_test == IERSConventions.COEFFS_NUT1980_ψ
        @test NUT1980_ϵ_test == IERSConventions.COEFFS_NUT1980_ϵ

    end

end;
