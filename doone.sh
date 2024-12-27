avid=$1
p=$2
printf "%s %si | ./doone avid=$avid -p=$p launched.\n" "$(date)" "$line"
printf "%s %si | ./doone avid=$avid -p=$p launched.\n" "$(date)" "$line" >> log.log

cd AVAC

make run "avid=$avid"
ec=$?
printf "%s %si | TSUL/_output$avid exit code: $ec\n" "$(date)" "$line"
printf "%s %si | TSUL/_output$avid exit code: $ec\n" "$(date)" "$line" >> ../log.log
 
if [ "$p" == "-p" -a $? -eq 0 ]
then
    make iplot OUTDIR="_output$avid"
fi

cd -
cd TSUL # cd ../TSUL

if [ $ec -eq 0 ]
then
    make run "avid=$avid"
fi
ec=$?
printf "%s %si | AVAC/_output$avid exit code: $ec\n" "$(date)" "$line"
printf "%s %si | AVAC/_output$avid exit code: $ec\n" "$(date)" "$line" >> ../log.log

if [ "$p" == "-p" -a $ec -eq 0 ]
then
    make iplot OUTDIR="_output$avid"
fi

cd -
if [ "$p" == "-p" -a $ec -eq 0 ]
then
    python energy.py $avid
fi
tail log.log
