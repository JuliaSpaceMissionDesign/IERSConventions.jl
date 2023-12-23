
r2a = 180 / π * 3600

@testset "ERA" verbose=true begin 

    # Test on ERA angles (< 1 μas)
    for _ in 1:50

        tt_c = rand()/2
        tt_d = tt_c*Tempo.CENTURY2DAY

        ERA = iers_era(iers2010a, tt_d)
        ERAₑ = era00(DJ2000, tt_d)

        @test r2a*abs(ERA - ERAₑ) ≤ 1e-7
    end

end;