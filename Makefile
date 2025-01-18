install:
	sudo apt install libgdal-dev
	pip install -r pyrequirements.txt
	pip install --src=$(CLAW)/.. --no-build-isolation -e git+https://github.com/clawpack/clawpack.git@v5.11.0#egg=clawpack
	echo "PETSC_DIR/PETSC_ARCH/lib/... could match:"
	find / -name petsc 2>/dev/null

topo:
	$(CLAW_PYTHON) TOPM/topo.py

new:
	cd AVAC; make new
	cd TSUL; make new

avid ?= ""
one:
	cd AVAC; make run avid=$(avid)
	cd TSUL; make run avid=$(avid)

RUN_COUNT?=25
all:
	for avid in $(shell seq 4 $(RUN_COUNT)) ; do \
        make one avid=$$avid ; \
		python energy.py $$avid -s ; \
		git add -f log.log && git add figures/*.pdf && git commit -m "'make one avid=$(avid)' terminated" && git push ; \
    done

figures:
	for avid in $(shell seq 4 $(RUN_COUNT)) ; do \
        python energy.py -s $$avid ; \
    done

