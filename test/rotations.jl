r2a = 180 / π * 3600
v2as = (x, y) -> acosd(max(-1, min(1, dot(x / norm(x), y / norm(y))))) * 3600

@testset "Full Rotation without EOP" verbose=true begin 
    for _ in 1:100 

        tt_c = rand() 
        tt_d = tt_c*Tempo.CENTURY2DAY
        
        v = rand(BigFloat, 3)
        v /= norm(v)

        # TODO: write these tests 

    end
end;

@testset "Full Rotation with EOP" verbose=true begin 
    # TODO: write these tests 

end;
