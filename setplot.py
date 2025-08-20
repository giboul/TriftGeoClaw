#!/usr/bin/env python
""" 
Set up the plot figures, axes, and items to be done for each frame.

This module is imported by the plotting routines and then the
function setplot is called to set the plot parameters.
"""
from matplotlib import pyplot as plt
from pathlib import Path
import numpy as np
import utils
from clawpack.visclaw.data import ClawPlotData
from clawpack.visclaw import geoplot, gaugetools, plot_timing_stats


def mask_coarse(current_data):
    patch = current_data.framesoln.state.patch
    xc_centers,yc_centers = patch.grid.c_centers
    mask_coarse = np.empty(xc_centers.shape, dtype=bool)
    mask_coarse.fill(False)
    for state_fine in current_data.framesoln.states:
        # iterate over all patches, and find any finer level grids that are
        # sitting on top of this patch/grid/state.
        patch_fine = state_fine.patch
    
        # Only look at patches one level finer
        if patch_fine.level != patch.level+1:
            continue
    
        xlower_fine = patch_fine.dimensions[0].lower
        xupper_fine = patch_fine.dimensions[0].upper
        ylower_fine = patch_fine.dimensions[1].lower
        yupper_fine = patch_fine.dimensions[1].upper
    
        m1 = (xc_centers > xlower_fine) & (xc_centers < xupper_fine)
        m2 = (yc_centers > ylower_fine) & (yc_centers < yupper_fine)
    
        # Mask all fine grid regions
        mask_coarse = (m1 & m2) | mask_coarse
    
    current_data.add_attribute('mask_coarse',mask_coarse)


def setplot(plotdata: ClawPlotData = None) -> ClawPlotData:
    """ 
    Specify what is to be plotted at each frame.
    Input:  plotdata, an instance of pyclaw.plotters.data.ClawPlotData.
    Output: a modified version of plotdata.
    """ 

    utils.set_transparent_cmaps()
    background, back_extent = utils.read_world_image("topo_big.png")
    dam = utils.read_geojson("dam.geojson")[1:]
    def background_image(_):
        plt.imshow(background, extent=back_extent, zorder=0)
        plt.fill(*dam, c="k", zorder=0)

    if plotdata is None:
        plotdata = ClawPlotData()

    plotdata.clearfigures()  # clear any old figures,axes,items data
    plotdata.format = utils.config.get("output_format", "ascii")  # 'ascii' or 'binary' to match setrun.py

    # To plot gauge locations on pcolor or contour plot, use this as
    # an afteraxis function:
    def addgauges(current_data):
        gaugetools.plot_gauge_locations(current_data.plotdata,
                                        gaugenos='all',
                                        format_string='ko',
                                        add_labels=True)

    #-----------------------------------------
    # Figure for surface
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='Surface', figno=0)
    plotfigure.use_for_kml = True
    clawdata = utils.read_clawdata(Path(plotdata.outdir) / "claw.data")
    xmin, ymin = clawdata["lower"]
    xmax, ymax = clawdata["upper"]
    plotfigure.kml_xlimits = xmin, xmax
    plotfigure.kml_ylimits = ymin, ymax

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes('pcolor')
    plotaxes.title = 'Surface'
    plotaxes.scaled = True

    def fixup(current_data):

        t = current_data.t
        plt.gca().set_title(f'$hu$ at {t//60:.0f}min {t%60:.2f}s', fontsize=20)
        # plt.title("")
        # plt.xticks(fontsize=15)
        # plt.yticks(fontsize=15)

    plotaxes.beforeaxes = background_image
    plotaxes.afteraxes = fixup

    # Water
    def masked_var(data):
        # mask_coarse(data)
        drytol = data.user.get("dry_tolerance", 1e-3)
        h, hu, hv, eta = np.ma.masked_where(np.tile(data.q[0], (4,1,1))<=drytol, data.q)
        # arr = np.sqrt(hu**2 + hv**2)
        arr = eta - np.maximum(1767, eta-h)
        if hasattr(data, "mask_coarse"):
            return np.ma.masked_where(data.mask_coarse, arr)
        return arr

    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    # plotitem.plot_var = geoplot.surface_or_depth
    plotitem.plot_var = masked_var
    # plotitem.plot_var = geoplot.surface_or_depth
    plotitem.pcolor_cmap = "RdBu_water"
    plotitem.pcolor_cmin = -3 #topoconfig["lake_alt"]-2
    plotitem.pcolor_cmax = 3 #topoconfig["lake_alt"]+2
    plotitem.add_colorbar = True
    # plotitem.amr_celledges_show = [1,1,1]
    # plotitem.patchedges_show = 1

    if 0:
        # Land
        plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
        plotitem.plot_var = geoplot.land
        plotitem.pcolor_cmap = plt.cm.viridis
        plotitem.pcolor_cmin = utils.config.get("lake_alt",  0.) - 120
        plotitem.pcolor_cmax = utils.config.get("lake_alt",  0.) + 380
        plotitem.add_colorbar = False
        # plotitem.amr_celledges_show = [0,0,0]
        # plotitem.patchedges_show = 1

        # add contour lines of bathy if desired:
        plotitem = plotaxes.new_plotitem(plot_type='2d_contour')
        plotitem.show = False
        plotitem.plot_var = geoplot.topo
        plotitem.contour_levels = np.linspace(-3000,-3000,1)
        plotitem.amr_contour_colors = ['y']  # color on each level
        plotitem.kwargs = {'linestyles':'solid','linewidths':2}
        plotitem.amr_contour_show = [1,0,0]  
        plotitem.celledges_show = 0
        plotitem.patchedges_show = 0


    #-----------------------------------------
    # Figures for gauges
    #-----------------------------------------
    # plotfigure = plotdata.new_plotfigure(name='Surface at gauges', figno=300, type='each_gauge')
    # plotfigure.clf_each_gauge = True

    # # Set up for axes in this figure:
    # plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = xmin, xmax
    plotaxes.ylimits = ymin, ymax
    # plotaxes.title = 'Surface'

    # # Plot surface as blue curve:
    # plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    # plotitem.plot_var = 3
    # plotitem.plotstyle = 'b-'

    # # Plot topo as green curve:
    # plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    # plotitem.show = False

    # def gaugetopo(current_data):
    #     q = current_data.q
    #     h = q[0,:]
    #     eta = q[3,:]
    #     topo = eta - h
    #     return topo
        
    # plotitem.plot_var = gaugetopo
    # plotitem.plotstyle = 'g-'

    # -----------------------------------------
    # Plots of timing (CPU and wall time):

    def make_timing_plots(plotdata):

        timing_plotdir = Path(plotdata.plotdir) / '_timing_figures'
        # system(f'mkdir -p {timing_plotdir}')
        timing_plotdir.mkdir(exist_ok=True)
        # adjust units for plots based on problem:
        units = dict(comptime="seconds", simtime="hours", cell="millions")
        plot_timing_stats.make_plots(outdir=plotdata.outdir, 
                                     make_pngs=True,
                                     plotdir=timing_plotdir, 
                                     units=units)

    otherfigure = plotdata.new_otherfigure(
        name='timing plots',
        fname='_timing_figures/timing.html'
    )
    otherfigure.makefig = make_timing_plots

    #-----------------------------------------
    # Parameters used only when creating html and/or latex hardcopy
    # e.g., via pyclaw.plotters.frametools.printframes:

    plotdata.printfigs = True                # print figures
    plotdata.print_format = 'png'            # file format
    plotdata.print_framenos = 'all'          # list of frames to print
    plotdata.print_gaugenos = []             # list of gauges to print
    plotdata.print_fignos = 'all'            # list of figures to print
    plotdata.html = True                     # create html files of plots?
    plotdata.html_homelink = '../README.html'   # pointer for top of index
    plotdata.latex = True                    # create latex file of plots?
    plotdata.latex_figsperline = 2           # layout of plots
    plotdata.latex_framesperline = 1         # layout of plots
    plotdata.latex_makepdf = False           # also run pdflatex?
    plotdata.parallel = True                 # make multiple frame png's at once

    return plotdata


def parse_args():
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("setplot", default="setplot.py", type=str, nargs="?")
    parser.add_argument("outdir", default="./_output", type=str, nargs="?")
    return parser.parse_args()

def main():
    args = parse_args()
    from clawpack.visclaw.Iplotclaw import Iplotclaw as IPC
    ip = IPC(args.setplot, args.outdir)
    ip.plotloop()


if __name__ == "__main__":
    main()
