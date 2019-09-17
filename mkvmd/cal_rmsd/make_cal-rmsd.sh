#!/bin/bash
# -args <ref> <pdb>/0 <dcd> <outname>
# First frame of the trajectory will be the reference structure
# i.e. if pdb!==0, ref=<pdb>; else ref=<first frame of dcd>

PSF=initial-noW.psf
TRJ=md.trr

if [ ! -f $TRJ ]; then
    echo -e "$TRJ \nTrajectory not found!"
    exit 0
fi

SELREF=("protein and name CA and not resid 420 to 448 640 to 655 1151 to 1167" \
    "protein and name CA and resid 24 to 338" \
    "protein and name CA and resid 339 to 590 and not resid 420 to 448 " \
    "protein and name CA and resid 1 to 23 591 to 661 762 to 891 and not resid 640 to 655" \
    "protein and name CA and resid 662 to 761" \
    "protein and name CA and resid 892 to 1077 1255 to 1300" \
    "protein and name CA and resid 1078 to 1254 and not resid 1151 to 1167")
SELRMSD=("protein and name CA and not resid 420 to 448 640 to 655 1151 to 1167" \
    "protein and name CA and resid 24 to 338" \
    "protein and name CA and resid 339 to 590 and not resid 420 to 448 " \
    "protein and name CA and resid 1 to 23 591 to 661 762 to 891 and not resid 640 to 655" \
    "protein and name CA and resid 662 to 761" \
    "protein and name CA and resid 892 to 1077 1255 to 1300" \
    "protein and name CA and resid 1078 to 1254 and not resid 1151 to 1167")
OUTPUT_NAME=("all" "rec1" "rec2" "wed" "pi" "ruvc" "nuc")

echo "" > vm_cal-rmsd.tcl
cat >> vm_cal-rmsd.tcl << EOF

set mol [mol new $PSF waitfor all]
mol addfile $TRJ waitfor all
set num_frames [molinfo \$mol get numframes]
set sel_all [atomselect top all]

EOF

for ii in {0..6}
do
    
    cat >> vm_cal-rmsd.tcl << EOF
set sel_ref0 [atomselect top "${SELREF[$ii]}" frame 0]
set sel_ref [atomselect top "${SELREF[$ii]}"]
set sel_rmsd0 [atomselect top "${SELRMSD[$ii]}" frame 0]
set sel_rmsd [atomselect top "${SELRMSD[$ii]}"]
set outfile [open "rmsd_${OUTPUT_NAME[$ii]}.dat" "w"]
for {set i 0} {\$i<\$num_frames} {incr i} {
    \$sel_all frame \$i
    \$sel_ref frame \$i
    \$sel_rmsd frame \$i
    \$sel_all move [measure fit \$sel_ref \$sel_ref0]
    set rmsd [measure rmsd \$sel_rmsd \$sel_rmsd0]
    puts \$outfile "\$i \t \$rmsd"
}
close \$outfile


EOF

done

cat >> vm_cal-rmsd.tcl << EOF
quit
EOF

vmd -dispdev text -e vm_cal-rmsd.tcl 
