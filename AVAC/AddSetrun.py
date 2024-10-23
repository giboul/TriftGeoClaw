#======= Topo =======
lake_alt = 1767
dam_alt = lake_alt + 10
flood_seed = 2670350, 1171500

bounds = dict(
    xmin = 2666500.,
    xmax = 2673700.,
    ymin = 1169500.,
    ymax = 1174200.,
)

resolution = 10

#======= Computation =======
nx = 40
ny = 40

nsim = 9
tmax = 90

dt_init = 1.

cfl_desired = 0.5
nb_max_iter = 500

refinement = 3
refinement_area = 0
DryWetLimit = 0.0001

nodatavalue = 1  # ??

out_format = "binary"
