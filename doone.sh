cd TOPM
python topo.py
cd ../AVAC
make new
make data
avid = $1
make qinit "avid=$avid"
make output "OUTDIR=_output$avid"
 
cd ../TSUL
make new
make qinit
make data

make flows "avid=$avid"
make data "avid=$avid"
make run "avid=$avid"
