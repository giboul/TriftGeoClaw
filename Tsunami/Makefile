# Makefile for Clawpack code in this directory.
# This version only sets the local files and is frequently changed
# options, and then includes the standard makefile pointed to by CLAWMAKE.
CLAWMAKE = $(CLAW)/clawutil/src/Makefile.common

# See the above file for details and a list of make options, or type
# $ make .help

CLAW_PKG = geoclaw                  # Clawpack package to use
EXE = xgeoclaw                 # Executable to create
SETRUN_FILE = setrun.py        # File containing function to make data
OUTDIR = _output               # Directory for output
SETPLOT_FILE = setplot.py      # File containing function to set plots
PLOTDIR = _plots               # Directory for plots


# Environment variable FC should be set to fortran compiler, e.g. gfortran

# ---------------------------------
# package sources for this program:
# ---------------------------------

GEOLIB = $(CLAW)/geoclaw/src/2d/shallow
include $(GEOLIB)/Makefile.geoclaw

# ---------------------------------------
# package sources specifically to exclude
# (i.e. if a custom replacement source 
#  under a different name is provided)
# ---------------------------------------

EXCLUDE_MODULES = \

EXCLUDE_SOURCES = \

# ----------------------------------------
# List of custom sources for this program:
# ----------------------------------------


MODULES = \

SOURCES = \
  ./bc2amr.f90 \
  $(CLAW)/riemann/src/rpn2_geoclaw.f \
  $(CLAW)/riemann/src/rpt2_geoclaw.f \
  $(CLAW)/riemann/src/geoclaw_riemann_utils.f \

#-------------------------------------------------------------------
# Include Makefile containing standard definitions and make options:
include $(CLAWMAKE)

# Construct the topography data
.PHONY: topo all
topo:
	$(CLAW_PYTHON) maketopo.py

plotloop:
	$(CLAW_PYTHON) ../plotloop.py

movie:
	$(CLAW_PYTHON) ../movie.py
	xdg-open movie.gif

.movie:
	$(MAKE) .plots
	$(CLAW_PYTHON) ../movie.py
	xdg-open movie.gif

all: 
	$(MAKE) clobber
	$(MAKE) clean
	$(MAKE) new
	$(MAKE) data
	$(MAKE) topo
	$(MAKE) output
	$(MAKE) .plots
	$(MAKE) .htmls
	$(MAKE) movie
