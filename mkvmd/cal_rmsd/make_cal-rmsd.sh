#!/bin/bash
# -args <ref> <pdb>/0 <dcd> <outname>
# First frame of the trajectory will be the reference structure
# i.e. if pdb!==0, ref=<pdb>; else ref=<first frame of dcd>

PSF=../../pdb2gmx/prot.pdb
PDB=0
DCD=/home/PHARMACY/chan.773/storage/cpf1/run_noDNA/md-1-80-pf100ps-fit.xtc
SELREF=("protein and resid 940 to 990")
SELRMSD=("protein and resid 790 to 860")
OUTPUT_NAME=("r790-860")

for ii in {0..0}
do
    ref=${SELREF[$ii]}
    rmsd=${SELRMSD[$ii]}
    outname=${OUTPUT_NAME[$ii]}
    sed -e "s/SELREF/$ref/g" \
        -e "s/SELRMSD/$rmsd/g" \
        template-vm-cal-rmsd-tcl > vm_cal-rmsd.tcl
    vmd -dispdev text -e vm_cal-rmsd.tcl -args $PSF $PDB $DCD $outname
done
