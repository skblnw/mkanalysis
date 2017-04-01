nframe=$1
outfile=$2
rm -i $outfile

for (( frame=0; frame<$nframe; frame++ ))
do
    echo -ne "Doing Frame $frame/$nframe\r"
    time=`echo "scale=3; $frame*0.05" | bc`
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

echo -ne "\nDone\n"
