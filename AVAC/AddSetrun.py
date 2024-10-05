#======= Topo =======
lake_level = 1767
dam_level = lake_level + 10
flood_seed = 2670300, 1171500

xmin = 2668500.
xmax  = 2671000.
ymin = 1170100.
ymax  = 1171900.

resolution = 10

#======= Computation =======
nx = 40
ny = 40

nsim = 50
tmax = 80

dt_init = 1.

cfl_desired = 0.5
nb_max_iter = 500

refinement = 3
refinement_area = 0
DryWetLimit = 0.0001

nodatavalue = 1  # ??

out_format = "binary"
