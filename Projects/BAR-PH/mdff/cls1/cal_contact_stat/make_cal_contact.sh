#!/bin/bash
FILENAME=output/cls1-gs3-water

if true; then
#if false; then
    rm -f $FILENAME*I1*.dat $FILENAME*I1*.full

    egrep -v "^#" LINES-I1 > TMP
    sed -e 's/NUMLIST/1 2 3 4 5 6 7 8/g' template-cal-contact > cal_contact.tcl
    while read line; do
        vmd -dispdev text -e cal_contact.tcl -args $FILENAME $line </dev/null
        wait
    done < TMP
    rm TMP

    sed -e 's/NUMLIST/1 2 3 4 5 6 7 8/g' template-cal-contact-all > cal_contact-all.tcl
    vmd -dispdev text -e cal_contact-all.tcl -args $FILENAME contact-I1-all
fi

if true; then
#if false; then
    rm -f $FILENAME*I2*.dat $FILENAME*I2*.full

    egrep -v "^#" LINES-I2 > TMP
    sed -e 's/NUMLIST/1 2 5 6 7 8 3 4/g' template-cal-contact > cal_contact.tcl
    while read line; do
        vmd -dispdev text -e cal_contact.tcl -args $FILENAME $line </dev/null
        wait
    done < TMP
    rm TMP

    sed -e 's/NUMLIST/1 2 5 6 7 8 3 4/g' template-cal-contact-all > cal_contact-all.tcl
    vmd -dispdev text -e cal_contact-all.tcl -args $FILENAME contact-I2-all
fi

rm LOG_total_contact_number.log
wc -l $FILENAME*.full  >> LOG_total_contact_number.log
tail $FILENAME*I1*.dat -n 1  >> LOG_total_contact_number.log
tail $FILENAME*I2*.dat -n 1  >> LOG_total_contact_number.log
