INPUT_NAME=LOG_output_com_12.log
outfile=./contact.dat
> $outfile

for frame in {0..10}
do
    echo -ne "Frame $frame/8426\r"
    printf "%0.2f\t" $frame >> $outfile
    for resid in {251..361}
    do
        res=`sed -n '/^Frame '$frame'$/,/Frame/p' $INPUT_NAME | awk '{print $4}' | grep "^$resid$" | wc -l`

        printf "%i\t" $res >> $outfile
    done
    printf "\n" >> $outfile
done
echo -e "\nFinished. Finally. Congraz.\n"
