avid ?= 1
RUN_COUNT?=26

topo:
	python topm/maketopo.py

new:
	cd avac; make new
	cd tsul; make new

qinit:
	python topm/makeqinit_avac.py $(avid)
	python topm/makeqinit_tsul.py

run:
	cd avac; make data && make qinit avid=$(avid) && make output
	cd tsul; make data && make output

all:
	make new
	make qinit
	for avid in $(shell seq 1 $(RUN_COUNT)) ; do \
        make run avid=$$avid ; \
		mv avac/_output avac/_output$$avid
		mv tsul/_output tsul/_output$$avid
		# python tsul/energy.py $$avid -s && \
		# git add -f log.log figures/*.gif && \
		# git add figures/*.pdf && \
		# git commit -m "'make run avid=$$avid' terminated" && \
		# git push ; \
    done

figures:
	for avid in $(shell seq 1 $(RUN_COUNT)) ; do \
        python energy.py -s $$avid ; \
    done

