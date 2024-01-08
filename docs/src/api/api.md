
# [Public Documentation](@id iers_api)

## [IERS Models](@id iers_models)

This is a list of the supported IERS models and their approximations that can be used to 
selected the algorithm associated to a specific IERS convention. 

### 2010 Conventions

```@docs 
iers2010a 
iers2010b 

CPNc 
CPNd
```

### 2003 Conventions 

```@docs 
iers2003a 
iers2003b 
```

### 1996 Conventions 

```@docs 
iers1996 
```

## Rotations 

### CIO-based rotations
```@docs 
iers_rot3_gcrf_to_cirf 
iers_rot3_gcrf_to_itrf
iers_rot3_gcrf_to_tirf 
iers_rot3_itrf_to_cirf 
iers_rot3_itrf_to_tirf
```

### Equinox-based rotations 
```@docs
iers_rot3_gcrf_to_gtod
iers_rot3_gcrf_to_mod 
iers_rot3_gcrf_to_pef 
iers_rot3_gcrf_to_tod 
iers_rot3_itrf_to_mod 
iers_rot3_itrf_to_tod 
iers_rot3_itrf_to_gtod 
iers_rot3_itrf_to_pef 
```

## Bias, Precession and Nutation 

```@docs 
iers_bias 
iers_nutation 
iers_nutation_comp
iers_obliquity 
iers_precession
iers_pb
iers_npb
```

## Celestial Intermediate Pole 

```@docs 
iers_cip_motion
cip_xy
cip_xys
cip_vector 
```

## Earth Rotation and Sidereal Time 

```@docs 
iers_era 
iers_era_rotm 

iers_gmst 
iers_gast 
```

## Polar Motion 

```@docs 
iers_polar_motion
```

## [EOP Data](@id eop_data)

```@docs 
eop_filename
eop_load_data! 
eop_generate_from_csv 
eop_generate_from_txt
```