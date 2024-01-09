
PrecompileTools.@setup_workload begin

    PrecompileTools.@compile_workload begin

        # Here we manually turn the EOP initialisation to true so 
        # that we are able to call all the routines. It doesn't really 
        # matter what values are returned as long as the functions do not throw errors.
        
        IERS_EOP.init = true     
        
        t = 0.0

        # Just calling this functions should be enough to pre-compile everything 
        # else that downstreams from those.
        for m in (iers2010a, iers2010b, CPNc, CPNd, iers2003a, iers2003b, iers1996)
            
            # Precompile all rotations!
            iers_rot3_gcrf_to_mod(t, m)
            iers_rot3_gcrf_to_tod(t, m)
            iers_rot3_gcrf_to_gtod(t, m)
            iers_rot3_gcrf_to_pef(t, m) 
            iers_rot3_itrf_to_mod(t, m) 
            iers_rot3_itrf_to_tod(t, m) 
            iers_rot3_itrf_to_gtod(t, m) 
            iers_rot3_itrf_to_pef(t, m)
            iers_rot3_gcrf_to_cirf(t, m) 
            iers_rot3_gcrf_to_itrf(t, m)
            iers_rot3_gcrf_to_tirf(t, m) 
            iers_rot3_itrf_to_cirf(t, m) 
            iers_rot3_itrf_to_tirf(t, m)

            iers_bias(m, t)
            iers_precession(m, t)

            iers_nutation(m, t)
            iers_nutation(m, t, t, t)

            precession_angles_rot3(m, t)
            precession_angles_rot4(m, t)

        end

        # We make sure that EOP are unloaded
        eop_unload_data!()

    end

end