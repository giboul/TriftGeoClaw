#!/usr/bin/env python
"""
Module to set up run time parameters for Clawpack.

The values set in the function setrun are then written out to data files
that will be read in by the Fortran code.
"""
from pathlib import Path
from argparse import ArgumentParser
import numpy as np
from clawpack.clawutil.data import ClawData
from clawpack.clawutil.data import ClawRunData
from clawpack.geoclaw.fgout_tools import FGoutGrid
from clawpack.geoclaw.fgmax_tools import FGmaxGrid
import bc_inflows
from config import config


def setrun(claw_pkg='geoclaw',
           AVAC_outdir: str=None,
           outdir="_output",
           bouss=False) -> ClawRunData:
    """
    Define the parameters used for running Clawpack.

    Output
    ------
        ClawRunData
    """

    if AVAC_outdir is None:
        raise ValueError("Please specify the AVAC output `AVAC_outdir` directory")
    AVAC_outdir = Path(AVAC_outdir).expanduser().absolute()
    avacprobdata = ClawData()
    avacprobdata.read(Path(AVAC_outdir) / "setprob.data", force=True)
    avacfgrid = FGoutGrid(1, Path(AVAC_outdir))
    avacfgrid.read_fgout_grids_data()

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

    # Lower and upper edge of computational domain:
    clawdata.lower = config["lower"]
    clawdata.upper = config["upper"]

    # Number of grid cells: Coarsest grid
    cell_size = config["cell_size"]
    coarse_dx = cell_size * np.prod(config["amr_ratios"]["x"])
    coarse_dy = cell_size * np.prod(config["amr_ratios"]["y"])
    clawdata.num_cells[0] = int((clawdata.upper[0]-clawdata.lower[0]) / coarse_dx)
    clawdata.num_cells[1] = int((clawdata.upper[1]-clawdata.lower[1]) / coarse_dy)

    # Check that bounds are inside the fixed grid
    # Else it can cause discontinuities, causing large CFLs
    assert (
        avacfgrid.x1 < clawdata.lower[0] and # -2*cell_size[0] and 
        avacfgrid.x2 > clawdata.upper[0] and # +2*cell_size[0] and 
        avacfgrid.y1 < clawdata.lower[1] and # -2*cell_size[1] and 
        avacfgrid.y2 > clawdata.upper[1] # +2*cell_size[1]
    ), (
        "AVAC fgout must contain the current bounds (including ghost cells).\n"
        f"{avacfgrid.x1:.1f}, {clawdata.lower[0]:.1f}, {clawdata.upper[0]:.1f}, {avacfgrid.x2:.1f}\n"
        f"{avacfgrid.y1:.1f}, {clawdata.lower[1]:.1f}, {clawdata.upper[1]:.1f}, {avacfgrid.y2:.1f}\n"
    )

    # ---------------
    # Size of system:
    # ---------------
    # Number of equations in the system:
    clawdata.num_eqn = 5 if bouss else 3
    # Number of auxiliary variables in the aux array (initialized in setaux)
    clawdata.num_aux = 1
    # Index of aux array corresponding to capacity function, if there is one:
    clawdata.capa_index = 0  # TODO: Add error 'mcapa should be nonzero' should mention 'capa_index'

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
        clawdata.num_output_times = config["nout"]
        clawdata.tfinal = clawdata.t0 + config["duration"]
        clawdata.output_t0 = True  # output at initial (or restart) time?

    elif clawdata.output_style == 2:
        # Specify a list of output times.
        clawdata.output_times = [0.5, 1.0]

    elif clawdata.output_style == 3:
        # Output every iout timesteps with a total of ntot time steps:
        clawdata.output_step_interval = 1
        clawdata.total_steps = 3
        clawdata.output_t0 = True

    clawdata.output_format = config["output_format"]   # 'ascii' or 'binary'
    clawdata.output_q_components = "all"   # h, hu, hv, eta
    # clawdata.output_aux_components = 'none'  # eta=h+B is in q
    clawdata.output_aux_onlyonce = False    # output aux arrays each frame

    # ---------------------------------------------------
    # Verbosity of messages to screen during integration:
    # ---------------------------------------------------
    # The current t, dt, and cfl will be printed every time step
    # at AMR levels <= verbosity.  Set verbosity = 0 for no printing.
    # (E.g. verbosity == 2 means print only on levels 1 and 2.)
    clawdata.verbosity = config["verbosity"]

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
    if config["mode"] == "src":
        clawdata.bc_lower[0] = 'wall'
        clawdata.bc_upper[0] = 'wall'
        clawdata.bc_lower[1] = 'wall'
        clawdata.bc_upper[1] = 'wall'
    elif config["mode"] == "bc":
        clawdata.bc_lower[0] = 'user'
        clawdata.bc_upper[0] = 'user'
        clawdata.bc_lower[1] = 'user'
        clawdata.bc_upper[1] = 'user'
    elif config["mode"] == "debug":
        print(f"NO FLUX INTRODUCTION FROM AVAC ({config['mode'] = })")
        clawdata.bc_lower[0] = 'wall'
        clawdata.bc_upper[0] = 'wall'
        clawdata.bc_lower[1] = 'wall'
        clawdata.bc_upper[1] = 'wall'
    else:
        raise ValueError(f"Inflow mode '{config['mode']}' is not 'bc' or 'src'")

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
    amrdata.refinement_ratios_x = config["amr_ratios"]["x"]
    amrdata.refinement_ratios_y = config["amr_ratios"]["y"]
    amrdata.refinement_ratios_t = config["amr_ratios"]["t"]

    coarse_cell_size = np.subtract(clawdata.upper, clawdata.lower) / clawdata.num_cells
    fine_cell_size = coarse_cell_size / np.prod([amrdata.refinement_ratios_x, amrdata.refinement_ratios_y], axis=1)
    print(f"Coarsest cells size : {coarse_cell_size[0]:.2f}, {coarse_cell_size[1]:.2f}")
    print(f"Finest cells size : {fine_cell_size[0]:.2f}, {fine_cell_size[1]:.2f}")

    # max number of refinement levels:
    if config["max_refinement"] >= 0:
        amrdata.amr_levels_max = config["max_refinement"]
    else:
        amrdata.amr_levels_max = len(amrdata.refinement_ratios_x)+1

    # --------------
    # Time stepping:
    # --------------
    # if dt_variable==1: variable time steps used based on cfl_desired,
    # if dt_variable==0: fixed time steps dt = dt_initial will always be used.
    clawdata.dt_variable = True  # 1 is a truthy value => same

    # Initial time step for variable dt.
    # If dt_variable==0 then dt=dt_initial for all steps:
    clawdata.dt_initial = 1.

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
    # to specify regions of refinement append lines of the form
    #  [minlevel,maxlevel,t1,t2,x1,x2,y1,y2]
    rundata.regiondata.regions = [[
        amrdata.amr_levels_max, amrdata.amr_levels_max,
        clawdata.t0, clawdata.tfinal,
        clawdata.lower[0], clawdata.upper[0],
        clawdata.lower[1], clawdata.upper[1],
    ]]

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
    fgout.x1 = clawdata.lower[0]#  - 10
    fgout.x2 = clawdata.upper[0]#  + 10
    fgout.y1 = clawdata.lower[1]#  - 10
    fgout.y2 = clawdata.upper[1]#  + 10
    fgout.nx = int((fgout.x2 - fgout.x1)/1.)
    fgout.ny = int((fgout.y2 - fgout.y1)/1.)
    fgout.tstart = clawdata.t0
    fgout.tend = clawdata.tfinal
    fgout.nout = clawdata.num_output_times
    fgout_grids.append(fgout)

    # == fgmax grids ==
    rundata.fgmax_data.num_fgmax_val = 5  # Save depth and speed see https://www.clawpack.org/fgmax.html
    fgmax_grids = rundata.fgmax_data.fgmax_grids

    fgmax = FGmaxGrid()
    fgmax.point_style = 2  # uniform rectangular x-y grid
    fgmax.x1 = fgout.x1
    fgmax.x2 = fgout.x2
    fgmax.y1 = fgout.y1
    fgmax.y2 = fgout.y2
    fgmax.dx = (fgout.x2-fgout.x1) / fgout.nx
    fgmax.tstart_max = fgout.tstart      # when to start monitoring max values
    fgmax.tend_max   = fgout.tend       # when to stop monitoring max values
    fgmax.dt_check   = (fgout.tend-fgout.tstart) / fgout.nout     # target time (sec) increment between updating
    fgmax.min_level_check = amrdata.amr_levels_max    # which levels to monitor max on
    fgmax.arrival_tol = 1.e-2    # tolerance for flagging arrival
    fgmax.interp_method = 0      # 0 ==> pw const in cells, recommended
    fgmax_grids.append(fgmax)    # written to fgmax_grids.data

    # == setprob data ==
    # Interacting with AVAC
    probdata = rundata.new_UserData(name='probdata',fname='setprob.data')
    probdata.add_param('mode', config["mode"], 'The method for introucing the avalanche')
    probdata.add_param('damping', float(avacprobdata.rho)/rundata.geo_data.rho, 'rho_snow/rho_water')
    probdata.add_param('min_alt_avac', float(config['min_alt_avac']),  'Lake altitude')
    probdata.add_param('AVAC_outdir', str(AVAC_outdir), 'The directory containing the fixed grid output of AVAC.')
    probdata.add_param('bc_size', int(config["bc_size"]), 'Number of cells to interpolate from for each boundary.')
    probdata.add_param('input_format', avacfgrid.output_format, 'Number of cells to interpolate from for each boundary.')
    probdata.add_param('fgout_fgno', int(config["fgout_fgno"]), 'Fixed grid output id to read in AVAC.')

    if config["mode"] == "bc":
        extent = clawdata.lower[0], clawdata.upper[0], clawdata.lower[1], clawdata.upper[1]
        bc_inflows.write(AVAC_outdir, extent, outdir, config["bc_size"])

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
    geo_data.sea_level = 0. if Path(config["qinit"][0]).is_file() else config["lake_alt"]
    rundata.geo_data.speed_limit = 100.
    
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
    topo_data.topofiles = [[config["topo_type"], Path(config["bathymetry"]).expanduser()]]

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
    if Path(config["qinit"][0]).is_file():
        rundata.qinit_data.qinit_type = 4
        rundata.qinit_data.qinitfiles = [config["qinit"]]
    else:
        print(f"\n\tWARNING: No qinit file was specified ! ({config['qinit'][0] = })\n")

    if bouss is True:
        print("Adding BoussData")
        from clawpack.geoclaw.data import BoussData
        rundata.add_data(BoussData(), 'bouss_data')

        rundata.bouss_data.bouss_equations = 2    # 0=SWE, 1=MS, 2=SGN
        rundata.bouss_data.bouss_min_level = 1    # coarsest level to apply bouss
        rundata.bouss_data.bouss_max_level = 10   # finest level to apply bouss
        rundata.bouss_data.bouss_min_depth = 1.   # depth to switch to SWE
        rundata.bouss_data.bouss_solver = 3       # 1=GMRES, 2=Pardiso, 3=PETSc
        rundata.bouss_data.bouss_tstart = 0.      # time to switch from SWE

    return rundata


def parse_args():
    parser = ArgumentParser()
    parser.add_argument('claw_pkg', default='geoclaw', nargs='?')
    parser.add_argument('-a', '--AVAC_outdir', default=config["AVAC_outdir"])
    parser.add_argument('-b', '--bouss', action='store_true')
    parser.add_argument('-o', '--outdir', type=str, default="_output")
    return parser.parse_args()


def main():
    data = Path(".data")
    data.unlink(missing_ok=True)
    rundata = setrun(**parse_args().__dict__)
    rundata.write()
    data.touch()


if __name__ == '__main__':
    main()
