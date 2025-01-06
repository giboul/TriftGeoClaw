p=$1
printf "%s %si | TSUL/doall.sh -p=$p launched.\n" "$(date)" "$line"
printf "%s %si | TSUL/doall.sh -p=$p launched.\n" "$(date)" "$line" >> ../log.log

for ((avid = 0 ; avid < 24 ; avid++ ));
do
    make run "avid=$avid"
    if [ "$p" == "-p" -a $ec -eq 0 ]
    then
        make iplot OUTDIR="_output$avid"
        python energy.py $avid
    fi
    ec=$?
    printf "%s %si | TSUL/_output$avid exit code: $ec\n" "$(date)" "$line"
    printf "%s %si | TSUL/_output$avid exit code: $ec\n" "$(date)" "$line" >> ../log.log
done

tail ../log.log
