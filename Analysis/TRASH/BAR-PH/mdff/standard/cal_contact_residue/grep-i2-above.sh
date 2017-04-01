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
        resa=`awk '{print $2}' output/contact-i2a-$frame.dat | grep "^$resid$" | wc -l`
        resb=`awk '{print $2}' output/contact-i2b-$frame.dat | grep "^$resid$" | wc -l`
        res=$(expr $resa + $resa)

        printf "%i\t" $res >> $outfile
    done
    printf "\n" >> $outfile
done
