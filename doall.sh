#/bin/bash
for ((avid = 14 ; avid < 26 ; avid++ ));
do 
    if [ $? == 0 ]; then
        ./doone.sh $avid $1
    fi
done
