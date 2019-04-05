#!/bin/bash
FILENAME=output/cls1-6tetramer-gs3

if true; then
#if false; then
    rm -f $FILENAME*I1*.dat $FILENAME*I1*.full

    egrep -v "^#" LINES-I1 > TMP
    sed -e 's/NUMLIST/1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24/g' template-cal-contact > cal_contact.tcl
    while read line; do
        vmd -dispdev text -e cal_contact.tcl -args $FILENAME $line </dev/null
        wait
    done < TMP

    sed -e 's/NUMLIST/1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24/g' -e 's/and resid $seltext[1-2] //g' template-cal-contact > cal_contact-all.tcl
    vmd -dispdev text -e cal_contact-all.tcl -args $FILENAME contact-I1-all

    rm TMP cal_contact.tcl cal_contact-all.tcl
fi

if true; then
#if false; then
    rm -f $FILENAME*I2*.dat $FILENAME*I2*.full

    egrep -v "^#" LINES-I2 > TMP
    sed -e 's/NUMLIST/7 8 3 4 19 20 7 8 1 2 5 6 5 6 17 18 15 16 11 12 23 24 15 16 9 10 13 14 13 14 21 22/g' template-cal-contact > cal_contact.tcl
    while read line; do
        vmd -dispdev text -e cal_contact.tcl -args $FILENAME $line </dev/null
        wait
    done < TMP

    sed -e 's/NUMLIST/7 8 3 4 19 20 7 8 1 2 5 6 5 6 17 18 15 16 11 12 23 24 15 16 9 10 13 14 13 14 21 22/g' -e 's/and resid $seltext[1-2] //g' template-cal-contact > cal_contact-all.tcl
    vmd -dispdev text -e cal_contact-all.tcl -args $FILENAME contact-I2-all

    rm TMP cal_contact.tcl cal_contact-all.tcl
fi

rm LOG_contact_number.log
wc -l $FILENAME*.full >> LOG_contact_number.log
tail $FILENAME*I1*.dat -n 1 >> LOG_contact_number.log
tail $FILENAME*I2*.dat -n 1 >> LOG_contact_number.log
