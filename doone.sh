cd AVAC

avid=$1
make run "avid=$avid"
 
make iplot OUTDIR="_output$avid" && cd ../TSUL
make run "avid=$avid"
make iplot OUTDIR="_output$avid"
