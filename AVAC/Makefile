# Makefile for Clawpack code in this directory.
# This version only sets the local files and frequently changed
# options, and then includes the standard makefile pointed to by CLAWMAKE.
CLAWMAKE = $(CLAW)/clawutil/src/Makefile.common

# See the above file for details and a list of make options, or type
#   make .help
# at the unix prompt.


# Adjust these variables if desired:
# ----------------------------------

CLAW_PKG = geoclaw                  # Clawpack package to use
EXE = xgeoclaw                 # Executable to create
SETRUN_FILE = setrun.py        # File containing function to make data
OUTDIR = _output               # Directory for output
SETPLOT_FILE = setplot.py      # File containing function to set plots
PLOTDIR = _plots               # Directory for plots

FFLAGS = -O2 -fopenmp
OMP_NUM_THREADS = 15

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
  ./module_voellmy.f90 \
  ./qinit_module.f90 \

SOURCES = \
  ./src2.f90 \
  ./setprob.f90 \
  $(CLAW)/riemann/src/rpn2_geoclaw.f \
  $(CLAW)/riemann/src/rpt2_geoclaw.f \
  $(CLAW)/riemann/src/geoclaw_riemann_utils.f \

#-------------------------------------------------------------------
# Include Makefile containing standard definitions and make options:
include $(CLAWMAKE)

.PHONY: topo all
avid=""
data:
	-rm -f .data
	$(CLAW_PYTHON) $(SETRUN_FILE) $(CLAW_PKG) $(avid)
	touch .data
topo:
	$(CLAW_PYTHON) ../TOPM/topo.py
qinit:
	$(CLAW_PYTHON) qinit.py $(avid)
run:
	rm -f _output$(avid)/fgout*
	make data avid=$(avid) && make qinit avid=$(avid) && make output OUTDIR=_output$(avid)
	echo `date +'%y.%m.%d %H:%M:%S'` "| AVAC/_output$$avid exit code: $$?" >> ../log.log
flows:
	$(CLAW_PYTHON) flows.py $(avid)

RUN_COUNT?=26
runs:
	for avid in $(shell seq 1 $(RUN_COUNT)) ; do \
		make run avid=$$avid; \
	done
