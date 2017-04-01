outfile=./output2/contact-i1.dat
rm $outfile
for frame in {0..232}
do
    echo "Frame $frame"
    time=`echo "scale=2; $frame*0.1" | bc`
    printf "%0.2f\t" $time >> $outfile
    for resid in {1..361}
    do
        resa1=`awk '{print $1}' output/contact-i1a-$frame.dat | grep "^$resid$" | wc -l`
        resa2=`awk '{print $2}' output/contact-i1a-$frame.dat | grep "^$resid$" | wc -l`
        resb1=`awk '{print $1}' output/contact-i1b-$frame.dat | grep "^$resid$" | wc -l`
        resb2=`awk '{print $2}' output/contact-i1b-$frame.dat | grep "^$resid$" | wc -l`
        res=$(expr $resa1 + $resa2 + $resb1 + $resb2)

        printf "%i\t" $res >> $outfile
    done
    printf "\n" >> $outfile
done
