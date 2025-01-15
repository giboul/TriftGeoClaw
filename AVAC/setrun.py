#!/usr/bin/env python
"""
Module to set up run time parameters for Clawpack.

The values set in the function setrun are then written out to data files
that will be read in by the Fortran code.

"""
from argparse import ArgumentParser
from yaml import safe_load
from pathlib import Path
import numpy as np
from clawpack.geoclaw.fgout_tools import FGoutGrid


projdir = Path(__file__).parents[1]
with open(projdir / "config.yaml") as file:
    config = safe_load(file)
    TOPM = config["TOPM"]
    AVAC = config["AVAC"]
    TSUL = config["TSUL"]


#------------------------------
def setrun(claw_pkg='geoclaw', avid=""):
#------------------------------

    """
    Define the parameters used for running Clawpack.

    INPUT:
        claw_pkg expected to be "geoclaw" for this setrun.

    OUTPUT:
        rundata - object of class ClawRunData

    """

    from clawpack.clawutil import data

    assert claw_pkg.lower() == 'geoclaw',  "Expected claw_pkg = 'geoclaw'"

    num_dim = 2
    rundata = data.ClawRunData(claw_pkg, num_dim)

    #------------------------------------------------------------------
    # GeoClaw specific parameters:
    #------------------------------------------------------------------
    rundata = setgeo(rundata, avid)

    #------------------------------------------------------------------
    # Standard Clawpack parameters to be written to claw.data:
    #   (or to amr2ez.data for AMR)
    #------------------------------------------------------------------
    clawdata = rundata.clawdata  # initialized when rundata instantiated


    # ---------------
    # Spatial domain:
    # ---------------

    # Number of space dimensions:
    clawdata.num_dim = num_dim

    # Lower and upper edge of computational domain:
    bounds = TOPM["bounds"] | AVAC.get("bounds", dict())
    clawdata.lower[0] = bounds["xmin"]
    clawdata.upper[0] = bounds["xmax"]
    clawdata.lower[1] = bounds["ymin"]
    clawdata.upper[1] = bounds["ymax"]


    # Set single grid parameters first.
    # See below for AMR parameters.

    # Number of grid cells: Coarsest grid
    clawdata.num_cells[0] = int((clawdata.upper[0] - clawdata.lower[0])/66)
    clawdata.num_cells[1] = int((clawdata.upper[1] - clawdata.lower[1])/66)

    # ---------------
    # Size of system:
    # ---------------

    # Number of equations in the system:
    clawdata.num_eqn = 3

    # Number of auxiliary variables in the aux array (initialized in setaux)
    clawdata.num_aux = 1

    # Index of aux array corresponding to capacity function, if there is one:
    clawdata.capa_index = 0

    
    
    # -------------
    # Initial time:
    # -------------

    clawdata.t0 = 0.0


    # Restart from checkpoint file of a previous run?
    # If restarting, t0 above should be from original run, and the
    # restart_file 'fort.chkNNNNN' specified below should be in 
    # the OUTDIR indicated in Makefile.

    clawdata.restart = False               # True to restart from prior results
    clawdata.restart_file = 'fort.chk00006'  # File to use for restart data

    # -------------
    # Output times:
    #--------------

    # Specify at what times the results should be written to fort.q files.
    # Note that the time integration stops after the final output time.
    # The solution at initial time t0 is always written in addition.

    clawdata.output_style = 1

    if clawdata.output_style==1:
        # Output nout frames at equally spaced times up to tfinal:
        clawdata.num_output_times = 50
        clawdata.tfinal = 300
        clawdata.output_t0 = True  # output at initial (or restart) time?

    elif clawdata.output_style == 2:
        # Specify a list of output times.
        clawdata.output_times = [0.5, 1.0]

    elif clawdata.output_style == 3:
        # Output every iout timesteps with a total of ntot time steps:
        clawdata.output_step_interval = 1
        clawdata.total_steps = 1
        clawdata.output_t0 = True
        

    clawdata.output_format = AVAC['out_format']      # 'ascii' or 'binary' 

    clawdata.output_q_components = 'all'   # could be list such as [True,True]
    clawdata.output_aux_onlyonce = True    # output aux arrays only at t0


    # ---------------------------------------------------
    # Verbosity of messages to screen during integration:
    # ---------------------------------------------------

    # The current t, dt, and cfl will be printed every time step
    # at AMR levels <= verbosity.  Set verbosity = 0 for no printing.
    #   (E.g. verbosity == 2 means print only on levels 1 and 2.)
    clawdata.verbosity = 2



    # ------------------
    # Method to be used:
    # ------------------

    # Order of accuracy:  1 => Godunov,  2 => Lax-Wendroff plus limiters
    clawdata.order = 2
    
    # Use dimensional splitting? (not yet available for AMR)
    clawdata.dimensional_split = 'unsplit'
    
    # For unsplit method, transverse_waves can be 
    #  0 or 'none'      ==> donor cell (only normal solver used)
    #  1 or 'increment' ==> corner transport of waves
    #  2 or 'all'       ==> corner transport of 2nd order corrections too
    clawdata.transverse_waves = 2

    # Number of waves in the Riemann solution:
    clawdata.num_waves = 3
    
    # List of limiters to use for each wave family:  
    # Required:  len(limiter) == num_waves
    # Some options:
    #   0 or 'none'     ==> no limiter (Lax-Wendroff)
    #   1 or 'minmod'   ==> minmod
    #   2 or 'superbee' ==> superbee
    #   3 or 'mc'       ==> MC limiter
    #   4 or 'vanleer'  ==> van Leer
    # clawdata.limiter = ['mc', 'mc', 'mc']
    clawdata.limiter = ['minmod', 'minmod', 'minmod']

    clawdata.use_fwaves = True    # True ==> use f-wave version of algorithms
    
    # Source terms splitting:
    #   src_split == 0 or 'none'    ==> no source term (src routine never called)
    #   src_split == 1 or 'godunov' ==> Godunov (1st order) splitting used, 
    #   src_split == 2 or 'strang'  ==> Strang (2nd order) splitting used,  not recommended.
    clawdata.source_split = 'godunov'


    # --------------------
    # Boundary conditions:
    # --------------------

    # Number of ghost cells (usually 2)
    clawdata.num_ghost = 2

    # Choice of BCs at xlower and xupper:
    #   0 => user specified (must modify bcN.f to use this option)
    #   1 => extrapolation (non-reflecting outflow)
    #   2 => periodic (must specify this at both boundaries)
    #   3 => solid wall for systems where q(2) is normal velocity

    clawdata.bc_lower[0] = 'extrap'
    clawdata.bc_upper[0] = 'extrap'

    clawdata.bc_lower[1] = 'extrap'
    clawdata.bc_upper[1] = 'extrap'

    # Specify when checkpoint files should be created that can be
    # used to restart a computation.

    clawdata.checkpt_style = 0

    if clawdata.checkpt_style == 0:
        # Do not checkpoint at all
        pass

    elif abs(clawdata.checkpt_style) == 1:
        # Checkpoint only at tfinal.
        pass

    elif abs(clawdata.checkpt_style) == 2:
        # Specify a list of checkpoint times.  
        clawdata.checkpt_times = [0.1,0.15]

    elif abs(clawdata.checkpt_style) == 3:
        # Checkpoint every checkpt_interval timesteps (on Level 1)
        # and at the final time.
        clawdata.checkpt_interval = 5

    # ---------------
    # AMR parameters:
    # ---------------
    amrdata = rundata.amrdata
    amrdata.max1d = 60

    # List of refinement ratios at each level (length at least mxnest-1)
    amrdata.refinement_ratios_x = [2, 4, 4]
    amrdata.refinement_ratios_y = [2, 4, 4]
    amrdata.refinement_ratios_t = [2, 4, 4]

    # max number of refinement levels:
    amrdata.amr_levels_max = 1 + min(map(len, (amrdata.refinement_ratios_x, amrdata.refinement_ratios_y, amrdata.refinement_ratios_t)))

    cell_size_x = (clawdata.upper[0] - clawdata.lower[0])/clawdata.num_cells[0]
    min_cell_size_x = cell_size_x/np.prod(amrdata.refinement_ratios_x)
    print(f"Minimum cell size (x): {min_cell_size_x:.2f}")


    # --------------
    # Time stepping:
    # --------------

    # if dt_variable==1: variable time steps used based on cfl_desired,
    # if dt_variable==0: fixed time steps dt = dt_initial will always be used.
    clawdata.dt_variable = True

    # Initial time step for variable dt.
    # If dt_variable==0 then dt=dt_initial for all steps:
    clawdata.dt_initial = min_cell_size_x / 300.

    # Max time step to be allowed if variable dt used:
    clawdata.dt_max = 1e+99

    # Desired Courant number if variable dt used, and max to allow without
    # retaking step with a smaller dt:
    clawdata.cfl_desired = 0.5
    clawdata.cfl_max = 0.95

    # Maximum number of time steps to allow between output times:
    clawdata.steps_max = 500

    print(f"{clawdata.dt_initial = }")
    print(f"Max speed u: {min_cell_size_x/clawdata.dt_initial*clawdata.cfl_max:.2f} m/s")


    # Specify type of each aux variable in amrdata.auxtype.
    # This must be a list of length maux, each element of which is one of:
    #   'center',  'capacity', 'xleft', or 'yleft'  (see documentation).

    amrdata.aux_type = ['center']


    # Flag using refinement routine flag2refine rather than richardson error
    amrdata.flag_richardson = False    # use Richardson?
    amrdata.flag2refine = True

    # steps to take on each level L between regriddings of level L+1:
    amrdata.regrid_interval = 3

    # width of buffer zone around flagged points:
    # (typically the same as regrid_interval so waves don't escape):
    amrdata.regrid_buffer_width  = 2

    # clustering alg. cutoff for (# flagged pts) / (total # of cells refined)
    # (closer to 1.0 => more small grids may be needed to cover flagged cells)
    amrdata.clustering_cutoff = 0.700000

    # print info about each regridding up to this level:
    amrdata.verbosity_regrid = 0  


    #  ----- For developers ----- 
    # Toggle debugging print statements:
    amrdata.dprint = False      # print domain flags
    amrdata.eprint = False      # print err est flags
    amrdata.edebug = False      # even more err est flags
    amrdata.gprint = False      # grid bisection/clustering
    amrdata.nprint = False      # proper nesting output
    amrdata.pprint = False      # proj. of tagged points
    amrdata.rprint = False      # print regridding summary
    amrdata.sprint = False      # space/memory output
    amrdata.tprint = False      # time step reporting each level
    amrdata.uprint = False      # update/upbnd reporting
    
    # More AMR parameters can be set -- see the defaults in pyclaw/data.py

    # fgmax grid output for TSUL
    fgout_grids = rundata.fgout_data.fgout_grids  # empty list initially

    fgout = FGoutGrid()
    fgout.fgno = 1
    fgout.point_style = 2       # will specify a 2d grid of points
    xmin, xmax, ymin, ymax = expand_bounds(*np.loadtxt(projdir/"TOPM"/"lake_extent.txt"), 1/5)
    fgout.output_format = 'binary64'  # ascii, binary32 4-byte, float32
    # fgout.nx = int((xmax-xmin)/(clawdata.upper[0]-clawdata.lower[0]) * clawdata.num_cells[0]*np.prod(amrdata.refinement_ratios_x))
    # fgout.ny = int((ymax-ymin)/(clawdata.upper[1]-clawdata.lower[1]) * clawdata.num_cells[1]*np.prod(amrdata.refinement_ratios_y))
    fgout.nx = int((xmax-xmin)/1.)
    fgout.ny = int((ymax-ymin)/1.)
    fgout.x1 = xmin
    fgout.x2 = xmax
    fgout.y1 = ymin
    fgout.y2 = ymax
    fgout.tstart = 0.
    fgout.tend = clawdata.tfinal
    fgout.nout = clawdata.num_output_times
    fgout_grids.append(fgout)    # written to fgout_grids.data

    # == setregions.data values ==
    regions = rundata.regiondata.regions
    # to specify regions of refinement append lines of the form
    #  [minlevel,maxlevel,t1,t2,x1,x2,y1,y2]
    #regions.append([1, 1, 0., 1.e10, 925720,927160., 6451140.,6452155.])
       
    #regions.append([2, 3, 3., 1.e10,   52., 72.,   52., 72.])
    #regions.append([2, 3, 3., 1.e10,   75., 95.,   -10.,  10.])
    #regions.append([2, 4, 3.4, 1.e10,   57., 68.,   57., 68.])
    #regions.append([2, 4, 3.4, 1.e10,   83., 92.,   -4.,  4.])
 

    # == setgauges.data values ==
    # for gauges append lines of the form  [gaugeno, x, y, t1, t2]
    # rundata.gaugedata.add_gauge()

    # gauges along x-axis:
    # gaugeno = 0
    # for r in np.linspace(-10, 10., 10):
    #     gaugeno = gaugeno+1
    #     x = r + .001  # shift a bit away from cell corners
    #     y = .001
     #    rundata.gaugedata.gauges.append([gaugeno, x, y, 0., 1e10])

    friction = np.loadtxt(projdir / AVAC["friction"], skiprows=1)
    if avid == "":
        mu, xi, u_ = friction[:, 1:].mean(axis=0)
    else:
        mu, xi, u_ = friction[(friction[:, 0]==int(avid)).argmax(), 1:]
    voellmydata = rundata.new_UserData(name='probdata',fname='voellmy.data')
    voellmydata.add_param("snow_density", AVAC["snow_density"], "")
    voellmydata.add_param("xi", xi, "Voellmy: geometrical resistance")
    voellmydata.add_param("mu", mu, "Voellmy: friction coefficient ~snow viscosity")
    voellmydata.add_param("u_", u_, "Velocity threshold")
    voellmydata.add_param("beta_slope", 300.0, "Threshold bed slope")
    voellmydata.add_param("coulomb", 0, "Wether to use the Coulomb model")

    return rundata
    # end of function setrun
    # ----------------------


#-------------------
def setgeo(rundata, avid):
#-------------------
    """
    Set GeoClaw specific runtime parameters.
    For documentation see ....
    """

    if hasattr(rundata, 'geo_data'):
        geo_data = rundata.geo_data
    else:
        raise AttributeError("*** Error, this rundata has no 'geo_data' attribute")

    # == Physics ==
    geo_data.gravity = 9.81
    geo_data.coordinate_system = 1
    geo_data.earth_radius = 6367.5e3

    # == correctif val d'isère
    # rundata.topo_data.topo_missing = AddSetrun.nodatavalue

    # == Forcing Options
    geo_data.coriolis_forcing = False

    # == Algorithm and Initial Conditions ==
    geo_data.sea_level = 0
    geo_data.dry_tolerance = 1e-3
    geo_data.friction_forcing = True
    geo_data.manning_coefficient = 0.025
    geo_data.friction_depth = 20.0

    # Refinement data
    refinement_data = rundata.refinement_data
    refinement_data.wave_tolerance = 1.e-2
    refinement_data.variable_dt_refinement_ratios = True

    # == settopo.data values ==
    topo_data = rundata.topo_data
    # for topography, append lines of the form
    #    [topotype, minlevel, maxlevel, t1, t2, fname]
    # topo_data.topofiles.append([2, 1, 3, 0., 1.e10, 'topo.asc'])
    topo_data.topofiles = [[2, projdir/AVAC['topo']]]

    # == setdtopo.data values ==
    # dtopo_data = rundata.dtopo_data
    # for moving topography, append lines of the form :   (<= 1 allowed for now!)
    #   [topotype, minlevel,maxlevel,fname]

    # == setqinit.data values ==
    rundata.qinit_data.qinit_type = 1
    # for qinit perturbations, append lines of the form: (<= 1 allowed for now!)
    #   [minlev, maxlev, fname]
    rundata.qinit_data.qinitfiles = [[projdir/"AVAC"/f"qinit{avid}.xyz"]]

    # == setfixedgrids.data values ==
    # fixedgrids = rundata.fixed_grid_data.fixedgrids
    # for fixed grids append lines of the form
    # [t1,t2,noutput,x1,x2,y1,y2,xpoints,ypoints,\
    #  ioutarrivaltimes,ioutsurfacemax]

    return rundata
    # end of function setgeo
    # ----------------------


def expand_bounds(x1, x2, y1, y2, rel_margin=1/50, abs_margin=0):
    dx = (x2 - x1) * rel_margin + abs_margin
    dy = (y2 - y1) * rel_margin + abs_margin
    xmin = x1 - dx
    xmax = x2 + dx
    ymin = y1 - dy
    ymax = y2 + dy
    return xmin, xmax, ymin, ymax


def main():
    # Set up run-time parameters and write all data files.
    parser = ArgumentParser()
    parser.add_argument('claw_pkg', default='geoclaw', nargs='?')
    parser.add_argument('avid', default='', nargs='?')
    args = parser.parse_args()

    data = Path(".data")
    data.unlink(missing_ok=True)
    rundata = setrun(**args.__dict__)
    rundata.write(projdir / "AVAC")
    data.touch()


if __name__ == '__main__':
    main()

