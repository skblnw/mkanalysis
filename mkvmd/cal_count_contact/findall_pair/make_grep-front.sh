outfile=./contact_front.dat
rm -f $outfile
for frame in {0..81}
do
    echo -ne "Frame $frame/81\r"
    printf "%0.2f\t" $frame >> $outfile
    for resid in {1..361}
    do
        res1=`sed -n '/Frame '$frame'/,/Frame/p' LOG_findall-2.log | awk '{print $4}' | grep "^$resid$" | wc -l`
        res2=`sed -n '/Frame '$frame'/,/Frame/p' LOG_findall-2.log | awk '{print $8}' | grep "^$resid$" | wc -l`
        res=$(expr $res1 + $res2)

        printf "%i\t" $res >> $outfile
    done
    printf "\n" >> $outfile
done
echo -e "\nFinished. Finally. Congraz.\n"
