#!/bin/bash
# KevC @ 2019
# Bash script to calculate RMSD profile

PDB="$1"
TRJ="$2"
REF="$3"
SELREF=("chain B and name CA")
SELRMSD=("chain A and name CA")
[ $# -eq 0 ] && { echo "mkvmd> Usage: $0 [PDB] [TRJ] [REF]"; echo "mkvmd> By default, the selection is '$SELRMSD'"; exit 1; }

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

EOF

for ii in {0..0}
do
    jj=$(($ii+1))
    cat >> tcl << EOF
set sel_ref0 [atomselect 1 "${SELREF[$ii]}" frame 0]
set sel_ref [atomselect 0 "${SELREF[$ii]}"]
set sel_rmsd0 [atomselect 1 "${SELRMSD[$ii]}" frame 0]
set sel_rmsd [atomselect 0 "${SELRMSD[$ii]}"]
set outfile [open "$jj.dat" "w"]
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

cat >> tcl << 'EOF'
quit
EOF

vmd -dispdev text -e tcl

rm -f tcl