#======= Topo =======
scale = 1e3
sea_level = 1767/scale
# For setrun.py
cxstart = 2668850./scale
cxstop  = 2670894./scale
cystart = 1170285./scale
cystop  = 1171850./scale

xstart = 2668500./scale
xstop  = 2671000./scale
ystart = 1170100./scale
ystop  = 1171900./scale

if (
    cxstart <= xstart or
    cxstop >= xstop or
    cystart <= ystart or
    cystop >= ystop
):
    raise ValueError(
        "The cropping box (topo file) is not fully "
        "contained by the computation box (grid bounds)"
    )
# Normalize domain for computations
cxstart = cxstart - xstart
cxstop  = cxstop  - xstart 
cystart = cystart - ystart
cystop  = cystop  - ystart

xmin = cxstart
xmax = cxstop
ymin = cystart
ymax = cystop

#======= Computation =======
nx = 50
ny = 50

nsim = 50
tmax = 10

dt_init = 1

cfl_desired = 0.5
nb_max_iter = 500

refinement = 3
refinement_area = 0
DryWetLimit = 0.0001

nodatavalue = 1  # ??
