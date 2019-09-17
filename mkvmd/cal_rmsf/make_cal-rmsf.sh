#!/bin/bash
# -args <ref> <pdb>/0 <dcd> <outname>
# First frame of the trajectory will be the reference structure
# i.e. if pdb!==0, ref=<pdb>; else ref=<first frame of dcd>

BLOCK_RMSF=false

PSF=initial-noW.psf
PDB=initial-noW.pdb
TRJ=md.trr

if [ ! -f $TRJ ]; then
    echo -e "$TRJ \nTrajectory not found!"
    exit 0
fi


OUTPUT_NAME=("out")
REFSEL=("protein and name CA")

echo "" > vm_cal-rmsf.tcl
cat >> vm_cal-rmsf.tcl << EOF

set mol [mol new $PSF waitfor all]
mol addfile $PDB waitfor all
mol addfile $TRJ waitfor all
set num_frames [molinfo \$mol get numframes]
set seltext "$REFSEL"
set a [atomselect top \$seltext]
set ref [atomselect top \$seltext frame 0]
for {set ii 0} {\$ii < \$num_frames} {incr ii} {
    \$a frame \$ii
    \$a move [measure fit \$a \$ref]
}

EOF



if $BLOCK_RMSF; then

### Block RMSF starts here ###

    total_ns=2400
    psperframe=1200
    block_ns=240
    block_frame=$(($block_ns*1000/$psperframe))
    block_num=$(($total_ns/$block_ns))

    for ii in $(seq 0 $(($block_num-1)))
    do
        start_frame=$(($block_frame*$ii+1))
        stop_frame=$(($start_frame+$block_frame-1))
        echo $start_frame 
        echo $stop_frame

        cat >> vm_cal-rmsf.tcl << EOF

set outfile [open "rmsf_${OUTPUT_NAME}_$ii.dat" "w"]
set seltext "protein and name CA"
set sel [atomselect top \$seltext]
set rmsf [measure rmsf \$sel first $start_frame last $stop_frame]
for {set i 0} {\$i <  [\$sel num] } {incr i} {
    puts \$outfile "[expr {\$i+1}] [lindex \$rmsf \$i]"
}
close \$outfile

EOF
    done

    cat >> vm_cal-rmsf.tcl << EOF
    quit
EOF



else

### Ordinary RMSF starts here ###

    start_frame=0

    cat >> vm_cal-rmsf.tcl << EOF

set outfile [open "rmsf_$OUTPUT_NAME.dat" "w"]
set seltext "protein and name CA"
set sel [atomselect top \$seltext]
set rmsf [measure rmsf \$sel first $start_frame last -1]
for {set i 0} {\$i <  [\$sel num] } {incr i} {
    puts \$outfile "[expr {\$i+1}] [lindex \$rmsf \$i]"
}
close \$outfile

quit
EOF

fi

vmd -dispdev text -e vm_cal-rmsf.tcl 
