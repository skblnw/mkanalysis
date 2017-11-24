INPUT_NAME=LOG_output_com_95100_12.log
outfile=./avg.dat
> $outfile

START=6500
let NFRAME=8427-$START
for resid in {251..361}
do
    echo -ne "Residue $resid/361\r"
    printf "%i\t" $resid >> $outfile
    count=`sed -n '/^Frame '$START'$/,$p' $INPUT_NAME | grep "FOUND" | awk '{print $4}' | grep "^$resid$" | wc -l`
    res=`echo "scale=2; $count/$NFRAME" | bc`
    printf "%.2f\t" $res >> $outfile
    printf "\n" >> $outfile
done
echo -e "\nFinished. Finally. Congraz.\n"
