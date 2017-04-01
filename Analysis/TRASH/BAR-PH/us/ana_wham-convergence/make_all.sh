for i in 4
do
    j=`echo "$i*2" | bc`
    sh make_gwhamfile.sh $i $j
    wait
    k=`echo "$i*10+70" | bc`
    sh run_gwham.sh 65-${k}ns
    wait
done
