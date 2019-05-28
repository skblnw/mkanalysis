#!/bin/bash
# -args <ref> <pdb>/0 <dcd> <outname>
# First frame of the trajectory will be the reference structure
# i.e. if pdb!==0, ref=<pdb>; else ref=<first frame of dcd>

PSF=prot.pdb
PDB=0
DCD=(".xtc" ".xtc")

start_ns=0
psperframe=100
start_frame=$(($start_ns*1000/$psperframe))

OUTPUT_NAME=("out")
REFSEL=("protein and name CA")

for ii in {0..1}
do
    sel1=${REFSEL[0]}
    sed -e "s/REFSEL/$sel1/g" template-vm-cal-rmsf-tcl > vm_cal-rmsf.tcl
    vmd -dispdev text -e vm_cal-rmsf.tcl -args $PSF $PDB ${DCD[$ii]} $start_frame ${OUTPUT_NAME}-$ii
done
