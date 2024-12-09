#/bin/bash
for ((avid = 0 ; avid < 24 ; avid++ ));
do 
    if [ $? == 0 ]; then
        ./doone.sh $avid $1
    fi
done
