export iers_earth_rot_rate

# Nominal Earth angular velocity  
const ωₑ = 7.292_115_146_706_979e-5

""" 
    iers_earth_rot_rate(LOD::Number=0)

Compute the true angular velocity of the Earth accounting for the Length of the Day, i.e., 
the instantaneous rate of change of UT1 with respect to a uniform time scale. 
"""
function iers_earth_rot_rate(LOD::Number=0)
    return ωₑ * (1 - LOD / 86400)
end
