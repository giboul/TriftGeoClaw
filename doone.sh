cd AVAC

avid=$1
make run "avid=$avid"
 
if [[ $2 == "-p" ]]
then
    make iplot OUTDIR="_output$avid" &
fi

cd ../TSUL
make run "avid=$avid"

if [[ $2 == "-p" ]]
then
    make iplot OUTDIR="_output$avid" &
fi

cd ..
wait