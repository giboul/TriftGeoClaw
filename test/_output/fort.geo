  
 --------------------------------------------
 Physics Parameters:
 -------------------
    gravity:   9.8100000000000005     
    density water:   1025.0000000000000     
    density air:   1.1499999999999999     
    ambient pressure:   101300.00000000000     
    earth_radius:   6367500.0000000000     
    coordinate_system:           1
    sea_level:   0.0000000000000000     
  
    coriolis_forcing: F
    theta_0:   0.0000000000000000     
    friction_forcing: T
    manning_coefficient:   2.5000000000000001E-002
    friction_depth:   1000000000.0000000     
  
    dry_tolerance:   1.0000000000000001E-005
  
 --------------------------------------------
 Refinement Control Parameters:
 ------------------------------
    wave_tolerance:   9.9999999999999995E-007
    speed_tolerance:   1000000000000.0000        1000000000000.0000        1000000000000.0000        1000000000000.0000        1000000000000.0000        1000000000000.0000     
    Variable dt Refinement Ratios: T
 
  
 --------------------------------------------
 SETDTOPO:
 -------------
    num dtopo files =            0
  
 --------------------------------------------
 SETTOPO:
 ---------
    mtopofiles =            1
    
    /home/axel/Desktop/PDM/TriftGeoClaw/test/bathymetry.xyz                                                                                               
   itopotype =            1
   mx =           50   x = (  -1.5000000000000000      ,   1.5000000000000000      )
   my =           50   y = (  -1.5000000000000000      ,   1.5000000000000000      )
   dx, dy (meters/degrees) =    6.1224489795918366E-002   6.1224489795918366E-002
  
   Ranking of topography files  finest to coarsest:            1
  
  
 --------------------------------------------
 SETQINIT:
 -------------
 /home/axel/Desktop/PDM/TriftGeoClaw/test/qinit.xyz                                                                                                    
   
 Reading qinit data from
 /home/axel/Desktop/PDM/TriftGeoClaw/test/qinit.xyz                                                                                                    
   
  
 --------------------------------------------
 Multilayer Parameters:
 ----------------------
    check_richardson: T
    richardson_tolerance:  0.94999999999999996     
    eigen_method:           4
    inundation_method:           2
    dry_tolerance:   1.0000000000000001E-005
