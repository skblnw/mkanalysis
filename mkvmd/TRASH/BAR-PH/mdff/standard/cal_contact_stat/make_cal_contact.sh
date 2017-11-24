# vmd -dispdev text -e cal_contact.tcl -args <psfname> <dcdname> <outprefix> <seltext1> <seltext2> <outname>
# vmd -dispdev text -e cal_contacall.tcl -args <psfname> <dcdname> <outprefix> <outname>
psfname=../combine_dcd/cls1-4dimer.psf
dcdname=../combine_dcd/cls1-gs3-water-protein-8ns-s50.dcd
outprefix=output/cls1-gs3-water

if true; then
    rm -f $outprefix*I1*.dat $outprefix*I1*.full

    egrep -v "^#" LINES-I1 > TMP
    sed -e 's/NUMLIST/1 2 3 4 5 6 7 8/g' template-cal-contact > cal_contact.tcl
    while read line; do
        vmd -dispdev text -e cal_contact.tcl -args $psfname $dcdname $outprefix $line </dev/null
        wait
    done < TMP
    rm TMP

    sed -e 's/NUMLIST/1 2 3 4 5 6 7 8/g' template-cal-contact-all > cal_contact-all.tcl
    vmd -dispdev text -e cal_contact-all.tcl -args $psfname $dcdname $outprefix contact-I1-all
fi

if true; then
    rm -f $outprefix*I2*.dat $outprefix*I2*.full

    egrep -v "^#" LINES-I2 > TMP
    sed -e 's/NUMLIST/1 2 5 6 7 8 3 4/g' template-cal-contact > cal_contact.tcl
    while read line; do
        vmd -dispdev text -e cal_contact.tcl -args $psfname $dcdname $outprefix $line </dev/null
        wait
    done < TMP
    rm TMP

    sed -e 's/NUMLIST/1 2 5 6 7 8 3 4/g' template-cal-contact-all > cal_contact-all.tcl
    vmd -dispdev text -e cal_contact-all.tcl -args $psfname $dcdname $outprefix contact-I2-all
fi

wc -l $outprefix*.full
tail $outprefix*I1*.dat -n 1
tail $outprefix*I2*.dat -n 1
