#!/bin/bash
# Kev @ 2019
# Bash script to calculate RMSD profile

PDB="$1"
TRJ="$2"
OUTPUT="$3"
REF="$4"
# if [ ! -n "$REF" ]; then
#     REF="$PDB"
# fi
SELREF=("segname PROA and name CA and resid 423 to 538 564 to 640 710 to 746 766 to 807")
SELRMSD=("{segname PROA and backbone and resid 423 to 538 564 to 640 710 to 746 766 to 807} or {segname PROP and backbone and resid 1042 to 1174 1204 to 1233 1256 to 1300}")
[ $# -eq 0 ] && { echo "mkvmd> Usage: $0 [PDB] [TRJ] [OUTPUT] [REF]"; echo "mkvmd> By default, the selection is '$SELRMSD'"; exit 1; }

if [ ! -f $PDB ]; then
    echo -e "$PDB \nStructure not found!"
    exit 0
fi

if [ ! -f $TRJ ]; then
    echo -e "$TRJ \nTrajectory not found!"
    exit 0
fi

if [ ! -f $REF ]; then
    echo -e "$REF \nStructure not found!"
    exit 0
fi


cat > tcl << EOF
mol new $PDB waitfor all
mol addfile $TRJ waitfor all
set num_frames [molinfo top get numframes]
set sel_all [atomselect top all]
mol new $REF waitfor all
set sel_ref0 [atomselect 1 "${SELREF}" frame 0]
set sel_ref [atomselect 0 "${SELREF}"]
set sel_rmsd0 [atomselect 1 "${SELRMSD}" frame 0]
set sel_rmsd [atomselect 0 "${SELRMSD}"]
set outfile [open "$OUTPUT" "w"]
for {set i 0} {\$i<\$num_frames} {incr i} {
    \$sel_all frame \$i
    \$sel_ref frame \$i
    \$sel_rmsd frame \$i
    \$sel_all move [measure fit \$sel_ref \$sel_ref0]
    set rmsd [measure rmsd \$sel_rmsd \$sel_rmsd0]
    puts \$outfile "\$i \t \$rmsd"
}
close \$outfile
quit
EOF

vmd -dispdev text -e tcl
rm tcl
