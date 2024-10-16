cd AVAC

for ((avid = 0 ; avid < 24 ; avid++ ));
do 
    make qinit "avid=$avid"
    make output "OUTDIR=_output$avid"
    make flows "avid=$avid"
done

cd -
cd Tsunami
read -p "Do you wish to install this program? " yn
echo "Current working Directory: $(pwd)"
for ((avid = 0 ; avid < 24 ; avid++ ));
do
    make data "avid=$avid"
    make run "avid=$avid"
done
cd -
