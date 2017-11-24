#!/bin/bash
# -args <ref> <pdb>/0 <dcd> <outname>
# First frame of the trajectory will be the reference structure
# i.e. if pdb!==0, ref=<pdb>; else ref=<first frame of dcd>

PSF=../combine_dcd/initial-noW.psf
PDB=../combine_dcd/initial-noW.pdb
DCD=../combine_dcd/npt-pf1000ps.dcd

start_ns=100
psperframe=100
start_frame=$(($start_ns*1000/$psperframe))

OUTPUT_NAME=("P1-PROT" "P1-MONO" "P1-BAR" "P1-PH")
REFSEL=("protein and backbone" "protein and segname P1 and name CA" "protein and segname P1 and resid 1 to 250 and name CA" "protein and segname P1 and resid 251 to 361 and name CA")

for ii in {0..3}
do
    outname=${OUTPUT_NAME[$ii]}
    sel1=${REFSEL[$ii]}
    sel2=${RMSFSEL[$ii]}
    sed -e "s/REFSEL/$sel1/g" template-vm-cal-rmsf-tcl > vm_cal-rmsf.tcl
    vmd -dispdev text -e vm_cal-rmsf.tcl -args $PSF $PDB $DCD $start_frame $outname
done
