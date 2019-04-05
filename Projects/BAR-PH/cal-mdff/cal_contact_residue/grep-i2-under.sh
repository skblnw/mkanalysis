outfile=./output2/contact-i2-under.dat
rm -f $outfile
nframe=900

for frame in $(seq 1 $nframe)
do
    echo -ne "Frame $frame/$nframe\r"
    time=`echo "scale=2; $frame*0.1" | bc`
    printf "%0.2f\t" $time >> $outfile
    for resid in {1..361}
    do
        resa=`awk '{print $1}' output/contact-i2a-$frame.dat | grep "^$resid$" | wc -l`
        resb=`awk '{print $1}' output/contact-i2b-$frame.dat | grep "^$resid$" | wc -l`
        res=$(expr $resa + $resa)

        printf "%i\t" $res >> $outfile
    done
    printf "\n" >> $outfile
done
echo -e "\nFinished. Finally. Congraz.\n"
