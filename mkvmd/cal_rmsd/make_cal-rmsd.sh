#!/bin/bash
# -args <ref> <pdb>/0 <dcd> <outname>
# First frame of the trajectory will be the reference structure
# i.e. if pdb!==0, ref=<pdb>; else ref=<first frame of dcd>

PSF=../combine_dcd/initial-noW.psf
PDB=../combine_dcd/initial-noW.pdb
DCD=../combine_dcd/npt-pf1000ps.dcd
SELTEXT=("protein and backbone")
OUTPUT_NAME=("npt")

for ii in {0..0}
do
    tmp1=${SELTEXT[$ii]}
    tmp2=${OUTPUT_NAME[$ii]}
    sed -e "s/SELTEXT/$tmp1/g" template-vm-cal-rmsd-tcl > vm_cal-rmsd.tcl
    vmd -dispdev text -e vm_cal-rmsd.tcl -args $PSF $PDB $DCD $tmp2
done
