avid ?= ""

one:
	cd AVAC
	make run avid=$(avid)
	cd -
	cd ../TSUL
	make run avid=$(avid)
	cd -

all:
	for (( i=1; i<=23; i++ )) \
    do \
        make one avid=$$i ; \
    done
