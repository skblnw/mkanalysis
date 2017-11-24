for seg in {3..13}
do
    echo -e "\n\nComputing segment A$seg"
    outfile=./output2/contact-A$seg.dat
    rm $outfile
    for frame in {0..489}
    do
        echo -ne "Frame $frame/489\r"
        time=`echo "scale=2; $frame*0.5" | bc`
        printf "%0.2f\t" $time >> $outfile

        seg2=$(expr $seg - 2)
        res=`awk '{print $4}' output/contact-A$seg-$frame.dat | grep "^A${seg2}$" | wc -l`
        printf "%i\t" $res >> $outfile

        printf "\n" >> $outfile
    done
done
echo -e "\nFinished. Finally. Congraz.\n"
