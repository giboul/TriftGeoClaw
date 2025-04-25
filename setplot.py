#!/usr/bin/env python
""" 
Set up the plot figures, axes, and items to be done for each frame.

This module is imported by the plotting routines and then the
function setplot is called to set the plot parameters.
"""
from matplotlib import pyplot as plt
from pathlib import Path
from yaml import safe_load
import numpy as np

from clawpack.visclaw.data import ClawPlotData
from clawpack.visclaw import geoplot, gaugetools, plot_timing_stats

with open("config.yaml") as file:
    config = safe_load(file)

def setplot(plotdata: ClawPlotData = None) -> ClawPlotData:
    """ 
    Specify what is to be plotted at each frame.
    Input:  plotdata, an instance of pyclaw.plotters.data.ClawPlotData.
    Output: a modified version of plotdata.
    """ 

    if plotdata is None:
        plotdata = ClawPlotData()

    plotdata.clearfigures()  # clear any old figures,axes,items data
    plotdata.format = config["output_format"]  # 'ascii' or 'binary' to match setrun.py

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

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes('pcolor')
    plotaxes.title = 'Surface'
    plotaxes.scaled = True

    def fixup(current_data):

        # addgauges(current_data)
        t = current_data.t
        plt.gca().set_title(f'$hu$ at {t//60:.0f}min {t%60:.2f}s', fontsize=20)
        # plt.title("")
        # plt.xticks(fontsize=15)
        # plt.yticks(fontsize=15)

    plotaxes.afteraxes = fixup

    # Water
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = geoplot.surface
    # plotitem.plot_var = geoplot.surface_or_depth
    plotitem.plot_var = 1
    # plotitem.plot_var = geoplot.surface_or_depth
    plotitem.pcolor_cmap = plt.cm.Blues
    plotitem.pcolor_cmin = -10 #topoconfig["lake_alt"]-2
    plotitem.pcolor_cmax = 10 #topoconfig["lake_alt"]+2
    plotitem.add_colorbar = True
    # plotitem.amr_celledges_show = [1,1,1]
    # plotitem.patchedges_show = 1

    # Land
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = geoplot.land
    plotitem.pcolor_cmap = plt.cm.viridis
    plotitem.pcolor_cmin = config["lake_alt"] - 120
    plotitem.pcolor_cmax = config["lake_alt"] + 380
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
    # plotaxes.xlimits = 'auto'
    # plotaxes.ylimits = 'auto'
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


def main():
    from clawpack.visclaw.Iplotclaw import Iplotclaw as IPC
    ip = IPC()
    ip.plotloop()


if __name__ == "__main__":
    main()
