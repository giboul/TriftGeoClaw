#======= Topo =======
sea_level = 1767

bounds = dict(
    xmin = 2668500.,
    xmax  = 2671000.,
    ymin = 1170100.,
    ymax  = 1171900.,
)

resolution = 10
lake_alt = 1767
dam_alt = lake_alt + 10
flood_seed = (2.67026e6, 1.171490e6)
#======= Computation =======
nx = 40
ny = 40

nsim = 50
tmax = 70

dt_init = 1.

cfl_desired = 0.5
nb_max_iter = 500

refinement = 3
refinement_area = 0
DryWetLimit = 0.0001

nodatavalue = 1  # ??
out_format = "binary"
