#/bin/bash
cd TOPM && python topo.py && cd ../AVAC && make new && make data
for ((avid = 0 ; avid < 24 ; avid++ ));
do 
    if [ $? == 0 ]; then
        make qinit "avid=$avid" && make output "OUTDIR=_output$avid"
    fi
done

cd ../TSUL && make new && make qinit && make data

for ((avid = 0 ; avid < 24 ; avid++ ));
do
    if [ $? -eq 0 ]; then
        make run "avid=$avid"
    fi
done
cd ..
