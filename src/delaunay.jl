
# NOTE: 
# =================

# The Delaunay parameters for the 1996 conventions follow the ones of the latest SOFA 
# release rather than those published by the IERS. As stated by the Orekit team:
#
#   The differences are nevertheless quite small (about 4.8e-11 radians, which is 
#   sub-millimeter level for low Earth orbits).

""" 
    DelaunayArgs

Container holding the Delaunay arguments of the Nutation Theory. 

### Fields 
- `M` -- Mean anomaly of the Moon, in radians. 
- `S` -- Mean anomaly of the Sun, in radians. 
- `F` -- Mean longitude of the Moon minus mean longitude of the ascending node, in radians. 
- `D` -- Mean elongation of the Moon from the Sun, in radians. 
- `Ω` -- Mean longitude of the Moon's ascending node, in radians. 
"""
struct DelaunayArgs{T}
    M::T
    S::T
    F::T
    D::T
    Ω::T
end

function DelaunayArgs(m::IERSModel, tdb_c::Number)

    DelaunayArgs(
        delaunay_anomaly_moon(m, tdb_c), 
        delaunay_anomaly_sun(m, tdb_c), 
        delaunay_longitude_diff(m, tdb_c), 
        delaunay_elongation_moon(m, tdb_c), 
        delaunay_longitude_node(m, tdb_c)
    )

end 


"""
    delaunay_anomaly_moon(m::IERSModel, tdb_c::Number)

Compute the mean anomaly of the Moon, in radians, given time `tdb_c` expressed in TDB Julian 
centuries since J2000. 
"""
function delaunay_anomaly_moon(::IERSAModels, tdb_c::Number)
    # This expression is valid for both IERS2003 and IERS2010
    mod2pi(
        arcsec2rad(
            @evalpoly(tdb_c, 485868.249036, 1717915923.2178, 31.8792, 0.051635, -0.00024470)
        )
    )
end

function delaunay_anomaly_moon(::IERSModel, tdb_c::Number)
    arcsec2rad(mod(485868.249036 + 1717915923.2178tdb_c, 1296000))
end

function delaunay_anomaly_moon(::IERS1996, tdb_c::Number)
    r = mod(1325tdb_c, 1)
    mod2pi(arcsec2rad(@evalpoly(tdb_c, 485866.733, 715922.633, 31.310, 0.064)) + 2π*r)
end


"""
    delaunay_anomaly_sun(m::IERSModel, tdb_c::Number)

Compute the mean anomaly of the Sun, in radians, given time `tdb_c` expressed in TDB Julian 
centuries since J2000.
"""
function delaunay_anomaly_sun(::IERSAModels, tdb_c::Number)
    mod2pi(
        arcsec2rad(
            @evalpoly(tdb_c, 1287104.793048, 129596581.0481, -0.5532, 0.000136, -0.00001149)
        )
    )
end

function delaunay_anomaly_sun(::IERSModel, tdb_c::Number)
    arcsec2rad(mod(1287104.79305 + 129596581.0481tdb_c, 1296000))
end

function delaunay_anomaly_sun(::IERS1996, tdb_c::Number)
    r = mod(99tdb_c, 1)
    mod2pi(arcsec2rad(@evalpoly(tdb_c, 1287099.804, 1292581.224, -0.577, -0.012)) + 2π*r)
end


"""
    delaunay_longitude_diff(m::IERSModel, tdb_c::Number)

Compute the difference between the longitude of the Moon and the longitude of the Moon's 
node, in radians, given time `tdb_c` expressed in TDB Julian centuries since J2000.
"""
function delaunay_longitude_diff(::IERSAModels, tdb_c::Number)
    mod2pi(
        arcsec2rad(
            @evalpoly(tdb_c, 335779.526232, 1739527262.8478, -12.7512, -0.001037, +0.00000417)
        )
    )
end

function delaunay_longitude_diff(::IERSModel, tdb_c::Number)
    arcsec2rad(mod(335779.526232 + 1739527262.8478tdb_c, 1296000))
end

function delaunay_longitude_diff(::IERS1996, tdb_c::Number)
    r = mod(1342tdb_c, 1)
    mod2pi(arcsec2rad(@evalpoly(tdb_c, 335778.877, 295263.137, -13.257, 0.011)) + 2π*r)
end

"""
    delaunay_elongation_moon(m::IERSModel, tdb_c::Number)

Compute the mean elongation of the Moon from the Sun, in radians, given time `tdb_c` expressed in 
TDB Julian centuries since J2000.
"""
function delaunay_elongation_moon(::IERSAModels, tdb_c::Number)
    mod2pi(
        arcsec2rad(
            @evalpoly(tdb_c, 1072260.703692, 1602961601.2090, -6.3706, 0.006593, -0.00003169)
        )
    )
end

function delaunay_elongation_moon(::IERSModel, tdb_c::Number)
    arcsec2rad(mod(1072260.70369 + 1602961601.2090tdb_c, 1296000))
end

function delaunay_elongation_moon(::IERS1996, tdb_c::Number)
    r = mod(1236tdb_c, 1)
    mod2pi(arcsec2rad(@evalpoly(tdb_c, 1072261.307, 1105601.328, -6.891, 0.019)) + 2π*r)
end

"""
    delaunay_longitude_node(m::IERSModel, tdb_c::Number)

Compute the longitude of the mean ascending node of the lunar orbit on th ecliptic, 
measured from the mean equinox of date, in radians, given time `tdb_c` expressed in TDB Julian 
centuries since J2000.
"""
function delaunay_longitude_node(::IERSAModels, tdb_c::Number)
    mod2pi(
        arcsec2rad(
            @evalpoly(tdb_c, 450160.398036, -6962890.5431, 7.4722, 0.007702, -0.00005939)
        )
    )
end

function delaunay_longitude_node(::IERSModel, tdb_c::Number)
    arcsec2rad(mod(450160.398036 - 6962890.5431tdb_c, 1296000))
end

function delaunay_longitude_node(::IERS1996, tdb_c::Number)
    r = mod(-5tdb_c, 1)
    mod2pi(arcsec2rad(@evalpoly(tdb_c, 450160.280, -482890.539, 7.455, 0.008)) + 2π*r)
end