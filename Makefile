install:
	sudo apt install libgdal-dev
	pip install -r pyrequirements.txt
	pip install --src=$(CLAW)/.. --no-build-isolation -e git+https://github.com/clawpack/clawpack.git@v5.11.0#egg=clawpack
	echo "PETSC_DIR/PETSC_ARCH/lib/... could match:"
	find / -name petsc 2>/dev/null

topo:
	python topm/topo.py

new:
	cd avac; make new
	cd tsul; make new

avid ?= ""
one:
	cd avac; make qinit avid=$(avid) && make output && make iplot fps=5 fname=../figures/avac$(avid).gif
	cd tsul; make data avid=$(avid) && make output && make iplot fps=5 fname=../figures/tsul$(avid).gif

RUN_COUNT?=26
all:
	for avid in $(shell seq 5 $(RUN_COUNT)) ; do \
		rm -rv **/_output*; \
        make one avid=$$avid ; \
		python tsul/energy.py $$avid -s ; \
		git add -f log.log figures/*.gif && git add figures/*.pdf && git commit -m "'make one avid=$$avid' terminated" && git push ; \
    done

figures:
	for avid in $(shell seq 5 $(RUN_COUNT)) ; do \
        python energy.py -s $$avid ; \
    done

