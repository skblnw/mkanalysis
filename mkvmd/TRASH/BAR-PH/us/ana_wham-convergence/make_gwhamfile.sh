#!/bin/bash

initstart=15
#initend="$2"
initdiff1="$1"
initdiff2="$2"
suffix=full
foldername=../../output


echo -n '' > pdofilelist.dat
for i in 616 630 640 650 660 670 680 695 710 725 740 755 765 770 785 800 815 830 845 860 875 890 905 920 935 950 965 980
do
    if [[ {616,630,640,650,660,670} =~ $i ]]; then
        startj=$initstart
        endj=`expr $startj + $initdiff1`
        echo $startj $endj
    elif [[ {680,695,710,725,740,755,765,770,785} =~ $i ]]; then
        startj=`expr $initstart - 2`
        endj=`expr $startj + $initdiff1`
        echo $startj $endj
    else
        startj=`expr $initstart - 2`
        endj=`expr $startj + $initdiff2`
        echo $startj $endj
    fi

    fname=us-z$i
    echo -n '' > $fname.$suffix
    for (( j=$startj; j<=$endj; j++ )); do
        cat $foldername/$fname-$j.colvars.traj >> $fname.$suffix
    done

    wincenter=`echo "scale=2; $i*0.01" | bc` # gromacs uses nm; scale=2 gives two digits after the floating point.
    egrep -v "^\#" $fname.$suffix | awk '{time=(NR*0.04); deviation=($2-'$i'*0.1)*0.1; print time, deviation}' > data # gromacs uses ps; requires deviation from the umbrella window center instead of actual z
    sed -e 's/POSnm/'$wincenter'/' template-header > header # generate the header required by g_wham, note that force constant is set to 4.0 kcal/mol/AA (1637.6 kJ/mol/nm/nm) in template-header
    cat header data > $fname.dat
    echo "$fname.dat" >> pdofilelist.dat # contains file names of all data generated above, to be fed to g_wham
    echo -e "Done for $fname"
done
rm data
rm header
rm *.$suffix
