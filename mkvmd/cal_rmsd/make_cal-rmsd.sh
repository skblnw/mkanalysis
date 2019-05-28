#!/bin/bash
# -args <ref> <pdb>/0 <dcd> <outname>
# First frame of the trajectory will be the reference structure
# i.e. if pdb!==0, ref=<pdb>; else ref=<first frame of dcd>

PSF=ionized.pdb
PDB=0
DCD=(".xtc" ".xtc")

SELREF=("protein and name CA")
SELRMSD=("protein and name CA")
OUTPUT_NAME=("out")

for ii in {0..1}
do
    ref=${SELREF[0]}
    rmsd=${SELRMSD[0]}
    sed -e "s/SELREF/$ref/g" \
        -e "s/SELRMSD/$rmsd/g" \
        template-vm-cal-rmsd-tcl > vm_cal-rmsd.tcl
    vmd -dispdev text -e vm_cal-rmsd.tcl -args $PSF $PDB ${DCD[$ii]} ${OUTPUT_NAME}-$ii
done
