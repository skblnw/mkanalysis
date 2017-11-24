#!/bin/bash

read -p "Output folder name? [output-XXX]" name
OUT_DIR=./output-$name
mkdir -p $OUT_DIR

if [ -s $OUT_DIR/LOG_findall.log ]; then
    read -p "$OUT_DIR/LOG_findall.log already exists. Find all contact again?" -n 1 -r
    echo 
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        vmd -e vm_measure-findall.tcl | tee $OUT_DIR/LOG_findall.log
        wait
        echo -e "Done finding all"
    fi
else
    vmd -e vm_measure-findall.tcl | tee $OUT_DIR/LOG_findall.log
    wait
    echo -e "Done finding all"
fi

grep FOUND $OUT_DIR/LOG_findall.log | awk '{print $4" "$8}' | sort -u > $OUT_DIR/PAIRS
pairs=`tr '\n' '\ ' < $OUT_DIR/PAIRS`
echo -e "Done getting a list of residue pairs: $pairs"
echo -e "Start measuring contacts for them"

sed -e 's/ATOMPAIRS/'"$pairs"'/g' template-vm-measure-pair-tcl > vm_measure-pair.tcl
vmd -e vm_measure-pair.tcl | tee $OUT_DIR/LOG_pair.log
wait
echo -e "Done. Check up $OUT_DIR/ !"
