#!/usr/bin/env python
"""
Module to set up run time parameters for Clawpack.

The values set in the function setrun are then written out to data files
that will be read in by the Fortran code.
"""
from subprocess import Popen
from argparse import ArgumentParser
from yaml import safe_load
from pathlib import Path
import numpy as np
from clawpack.clawutil.data import ClawRunData
from clawpack.geoclaw.fgout_tools import FGoutGrid


projdir = Path(__file__).parents[1]
with open(projdir / "config.yaml") as file:
    config = safe_load(file)
    AVAC = config["AVAC"]
    TOPM = config["TOPM"]
    TSUL = config["TSUL"]

def setrun(claw_pkg='geoclaw', AVAC_DIR: str="_output", bouss=False, inflow="bc") -> ClawRunData:
    """
    Define the parameters used for running Clawpack.

    Output
    ------
        ClawRunData
    """
    num_dim = 2
    rundata = ClawRunData(claw_pkg, num_dim)
    rundata = setgeo(rundata, bouss)
    # Standard Clawpack parameters to be written to claw.data:
    # (or to amr2ez.data for AMR)
    clawdata = rundata.clawdata  # initialized when rundata instantiated

    # Set single grid parameters first.
    # See below for AMR parameters.

    #----------------
    # Spatial domain:
    #----------------
    # Number of space dimensions:
    clawdata.num_dim = num_dim

    xmin, xmax, ymin, ymax = np.loadtxt(projdir/"topm"/"lake_extent.txt")
    # Number of grid cells: Coarsest grid
    clawdata.num_cells[0] = int((xmax-xmin)/60)
    clawdata.num_cells[1] = int((ymax-ymin)/60)

    # Lower and upper edge of computational domain:
    clawdata.lower[0] = float(xmin)
    clawdata.upper[0] = float(xmax)
    clawdata.lower[1] = float(ymin)
    clawdata.upper[1] = float(ymax)

    # ---------------
    # Size of system:
    # ---------------
    # Number of equations in the system:
    clawdata.num_eqn = 5 if bouss else 3
    # Number of auxiliary variables in the aux array (initialized in setaux)
    clawdata.num_aux = 1
    # Index of aux array corresponding to capacity function, if there is one:
    clawdata.capa_index = 0  # TODO: Add comment to error 'mcapa should be nonzero'

    # -------------
    # Initial time:
    # -------------
    clawdata.t0 = 0.
    # Restart from checkpoint file of a previous run?
    # If restarting, t0 above should be from original run, and the
    # restart_file 'fort.chkNNNNN' specified below should be in 
    # the OUTDIR indicated in Makefile.
    clawdata.restart = False              # True to restart from prior results
    clawdata.restart_file = 'fort.chk00096'  # File to use for restart data

    # -------------
    # Output times:
    # -------------
    # Specify at what times the results should be written to fort.q files.
    # Note that the time integration stops after the final output time.
    # The solution at initial time t0 is always written in addition.
    clawdata.output_style = 1

    if clawdata.output_style==1:
        # Output nout frames at equally spaced times up to tfinal:
        clawdata.num_output_times = 100
        clawdata.tfinal = clawdata.t0 + 300
        clawdata.output_t0 = True  # output at initial (or restart) time?

    elif clawdata.output_style == 2:
        # Specify a list of output times.
        clawdata.output_times = [0.5, 1.0]

    elif clawdata.output_style == 3:
        # Output every iout timesteps with a total of ntot time steps:
        clawdata.output_step_interval = 1
        clawdata.total_steps = 3
        clawdata.output_t0 = True

    clawdata.output_format = TSUL["out_format"]   # 'ascii' or 'binary' 
    clawdata.output_q_components = "all"   # h, hu, hv, eta
    # clawdata.output_aux_components = 'none'  # eta=h+B is in q
    clawdata.output_aux_onlyonce = False    # output aux arrays each frame

    # ---------------------------------------------------
    # Verbosity of messages to screen during integration:
    # ---------------------------------------------------
    # The current t, dt, and cfl will be printed every time step
    # at AMR levels <= verbosity.  Set verbosity = 0 for no printing.
    # (E.g. verbosity == 2 means print only on levels 1 and 2.)
    clawdata.verbosity = 3

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
    clawdata.limiter = ['mc', 'mc', 'mc']

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
    if TSUL["inflow"] == "src":
        clawdata.bc_lower[0] = 'wall'
        clawdata.bc_upper[0] = 'wall'
        clawdata.bc_lower[1] = 'wall'
        clawdata.bc_upper[1] = 'wall'
    elif TSUL["inflow"] == "bc":
        clawdata.bc_lower[0] = 'user'
        clawdata.bc_upper[0] = 'user'
        clawdata.bc_lower[1] = 'user'
        clawdata.bc_upper[1] = 'user'
    else:
        raise ValueError(f"inflow mode '{TSUL['inflow']}' is not 'bc' or 'src'")

    # --------------
    # Checkpointing:
    # --------------
    # Specify when checkpoint files should be created that can be
    # used to restart a computation.
    clawdata.checkpt_style = 0

    if np.abs(clawdata.checkpt_style) == 2:
        # Specify a list of checkpoint times.  
        pass
        # clawdata.checkpt_times = [0.1,0.15]

    elif np.abs(clawdata.checkpt_style) == 3:
        # Checkpoint every checkpt_interval timesteps (on Level 1)
        # and at the final time.
        clawdata.checkpt_interval = 5

    # ---------------
    # AMR parameters:
    # ---------------
    amrdata = rundata.amrdata
    # maximum size of patches in each direction (matters in parallel):
    amrdata.max1d = 60

    # List of refinement ratios at each level (length at least mxnest-1)
    amrdata.refinement_ratios_x = [2, 4, 5]
    amrdata.refinement_ratios_y = [2, 4, 5]
    amrdata.refinement_ratios_t = [2, 4, 5]

    cell_size_x = (clawdata.upper[0] - clawdata.lower[0])/clawdata.num_cells[0]
    min_cell_size_x = cell_size_x / np.prod(amrdata.refinement_ratios_x)
    print(f"Minimum cells size (x): {min_cell_size_x:.2f}")

    # max number of refinement levels:
    amrdata.amr_levels_max = 4

    # --------------
    # Time stepping:
    # --------------
    # if dt_variable==1: variable time steps used based on cfl_desired,
    # if dt_variable==0: fixed time steps dt = dt_initial will always be used.
    clawdata.dt_variable = True  # 1 is a truthy value => same

    # Initial time step for variable dt.
    # If dt_variable==0 then dt=dt_initial for all steps:
    clawdata.dt_initial = float(min_cell_size_x) / 300.

    # Max time step to be allowed if variable dt used:
    clawdata.dt_max = 1e+99

    # Desired Courant number if variable dt used, and max to allow without
    # retaking step with a smaller dt:
    clawdata.cfl_desired = 0.75
    clawdata.cfl_max = 1.0

    # Maximum number of time steps to allow between output times:
    clawdata.steps_max = 5000

    # Specify type of each aux variable in amrdata.auxtype.
    # This must be a list of length maux, each element of which is one of:
    #   'center',  'capacity', 'xleft', or 'yleft'  (see documentation).
    amrdata.aux_type = ['center','capacity','yleft']

    # Flag using refinement routine flag2refine rather than richardson error
    amrdata.flag_richardson = False    # use Richardson?
    amrdata.flag_richardson_tol = 0.002  # Richardson tolerance
    amrdata.flag2refine = True

    # steps to take on each level L between regriddings of level L+1:
    amrdata.regrid_interval = 3

    # width of buffer zone around flagged points:
    # (typically the same as regrid_interval so waves don't escape):
    amrdata.regrid_buffer_width  = 2

    # clustering alg. cutoff for (# flagged pts) / (total # of cells refined)
    # (closer to 1.0 => more small grids may be needed to cover flagged cells)
    amrdata.clustering_cutoff = 0.7

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
    amrdata.tprint = True       # time step reporting each level
    amrdata.uprint = False      # update/upbnd reporting
    
    # More AMR parameters can be set -- see the defaults in pyclaw/data.py

    # --------
    # Regions:
    # --------
    rundata.regiondata.regions = []
    # to specify regions of refinement append lines of the form
    #  [minlevel,maxlevel,t1,t2,x1,x2,y1,y2]
    rundata.regiondata.regions.append([
        3, 4, clawdata.t0, clawdata.tfinal/5,
        2.67e6, 2.67e6+400,
        1.171e6+650, 1.171e6+950
    ])
    # -------
    # Gauges:
    # -------
    # rundata.gaugedata.gauges = []
    # for gauges append lines of the form  [gaugeno, x, y, t1, t2]
    # rundata.gaugedata.gauges.append([32412, xcoords.mean(), ycoords.mean(), clawdata.t0, tf])

    # == fgout grids ==
    # new style as of v5.9.0 (old rundata.fixed_grid_data is deprecated)
    # fixed_grid_data script doesn't exist anymore...

    fgout_grids = rundata.fgout_data.fgout_grids  # empty list initially

    fgout = FGoutGrid()
    fgout.fgno = 1
    fgout.point_style = 2       # will specify a 2d grid of points
    fgout.output_format = 'binary64'  # ascii, binary32 4-byte, float32
    fgout.nx = int((clawdata.upper[0]-clawdata.lower[0])/1.)
    fgout.ny = int((clawdata.upper[1]-clawdata.lower[1])/1.)
    fgout.x1 = clawdata.lower[0]
    fgout.x2 = clawdata.upper[0]
    fgout.y1 = clawdata.lower[1]
    fgout.y2 = clawdata.upper[1]
    fgout.tstart = clawdata.t0
    fgout.tend = clawdata.tfinal
    fgout.nout = clawdata.num_output_times
    fgout_grids.append(fgout)
 
    inflow_mode = TSUL.get("inflow", inflow)

    probdata = rundata.new_UserData(name='probdata',fname='setprob.data')
    probdata.add_param('mode', inflow_mode, 'The method for introucing the avalanche')
    probdata.add_param('damping', float(AVAC['snow_density'])/1025., 'rho_snow/rho_water')  # TODO water density
    probdata.add_param('lake_alt', float(TOPM['lake_alt']),  'Lake altitude')
    probdata.add_param('overhang', float(TOPM.get('overhang', 0.)), 'Overhang of the contour over the lake')
    if not Path(AVAC_DIR).is_absolute():
        AVAC_DIR = projdir / AVAC_DIR
    probdata.add_param('AVAC_DIR', str(AVAC_DIR), 'The directory containing the fixed grid output of AVAC.')

    return rundata


def setgeo(rundata: ClawRunData, bouss=False) -> ClawRunData:
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
    geo_data.sea_level = 0.

    # == Forcing Options ==
    geo_data.coriolis_forcing = False

    # == Algorithm and Initial Conditions ==
    geo_data.dry_tolerance = 1.e-10
    geo_data.friction_forcing = True
    geo_data.manning_coefficient =.025
    geo_data.friction_depth = 1e9

    # Refinement settings
    refinement_data = rundata.refinement_data
    refinement_data.variable_dt_refinement_ratios = True
    refinement_data.wave_tolerance = 1.e-6

    # == settopo.data values ==
    topo_data = rundata.topo_data
    # for topography, append lines of the form
    #    [topotype, fname]
    topo_data.topofiles = [[2, projdir / TOPM["bathymetry"]]]

    # == setdtopo.data values ==
    # dtopo_data = rundata.dtopo_data
    # # for moving topography, append lines of the form :   (<= 1 allowed for now!)
    # #   [topotype, fname]
    # dtopo_path = scratch_dir / 'dtopo_usgs100227.tt3'
    # dtopo_data.dtopofiles.append([3,dtopo_path])
    # dtopo_data.dt_max_dtopo = 0.2

    # == setqinit.data values ==
    # 0: nothing, 1: depth, 2: x momentum, 3: y momentum, 4: surface level
    # for qinit perturbations, append lines of the form: (<= 1 allowed for now!)
    #   [fname]
    # Check if using qinit or boundary condition
    rundata.qinit_data.qinit_type = 4
    rundata.qinit_data.qinitfiles = [["qinit.xyz"]]

    if bouss is True:
        print("Adding BoussData")
        from clawpack.geoclaw.data import BoussData
        rundata.add_data(BoussData(), 'bouss_data')
        
        rundata.bouss_data.bouss_equations = 2    # 0=SWE, 1=MS, 2=SGN # TODO
        rundata.bouss_data.bouss_min_level = 1    # coarsest level to apply bouss
        rundata.bouss_data.bouss_max_level = 10   # finest level to apply bouss
        rundata.bouss_data.bouss_min_depth = 1.   # depth to switch to SWE
        rundata.bouss_data.bouss_solver = 3       # 1=GMRES, 2=Pardiso, 3=PETSc
        rundata.bouss_data.bouss_tstart = 0.      # time to switch from SWE

    return rundata


def main():
    # Treating command line arguments
    parser = ArgumentParser()
    parser.add_argument('claw_pkg', default='geoclaw', nargs='?')
    parser.add_argument('AVAC_DIR', default='_output', nargs='?')
    parser.add_argument('--bouss', action='store_true')
    args = parser.parse_args()

    data = Path(".data")
    data.unlink(missing_ok=True)
    rundata = setrun(**args.__dict__)
    rundata.write(projdir / "tsul")
    data.touch()

    if TSUL["inflow"] == "bc":
        Popen(["python", projdir/"tsul"/"bc_inflows.py"])


def expand_bounds(x1, x2, y1, y2, rel_margin=1/50, abs_margin=0):
    dx = (x2 - x1) * rel_margin + abs_margin
    dy = (y2 - y1) * rel_margin + abs_margin
    xmin = x1 - dx
    xmax = x2 + dx
    ymin = y1 - dy
    ymax = y2 + dy
    return xmin, xmax, ymin, ymax

if __name__ == '__main__':
    main()
