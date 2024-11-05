#!/bin/bash
# KevC @ 2019
# Bash script to calculate RMSF profile 
# Turn on BLOCK_RMSF to calculate RMSF for blocks of trajectory

[ $# -eq 0 ] && { echo "mkvmd> Usage: $0 [PDB] [TRJ] [PREFIX] [START_FRAME] [REF]"; exit 1; }


PDB="$1"
TRJ="$2"
PREFIX="$3"
START_FRAME="$4"
REF="$5"

SELREF="noh chain A"
REGIONS=("name CA and chain A" "name CA and chain B" "name CA and chain C" "name CA and chain C and resid >303")
FILENAMES=("proa" "prob" "proc" "nhelix")

[ ! -z "PREFIX" ] && echo "mkvmd> Please specify prefix for output files"

if [ ! -n "$START_FRAME" ]; then
    START_FRAME=0
fi

if [ ! -n "$REF" ]; then
    REF="$PDB"
fi

files=("$PDB" "$TRJ" "$REF")
for file in "${files[@]}"; do
    if [ ! -f "$file" ]; then
        echo -e "$file \nStructure not found!"
        exit 1
    fi
done

cat > tcl <<EOF
mol new $PDB waitfor all
mol addfile $TRJ waitfor all
set num_frames [molinfo top get numframes]
set sel_all [atomselect top all]
mol new $REF waitfor all
set sel_ref0 [atomselect 1 "${SELREF}" frame 0]
set sel_ref [atomselect 0 "${SELREF}"]
EOF

counter=0
for region in "${REGIONS[@]}"; do
    filename="${PREFIX}-${FILENAMES[counter]}"
    rm $filename
    cat >> tcl <<EOF
set sel_rmsd0 [atomselect 1 "${region}" frame 0]
set sel_rms [atomselect 0 "${region}"]
set outfile [open $filename "w"]
puts \$outfile "# reference: $SELREF; calculated region: $region"
for {set i 0} {\$i<\$num_frames} {incr i} {
    \$sel_all frame \$i
    \$sel_ref frame \$i
    \$sel_rms frame \$i
    \$sel_all move [measure fit \$sel_ref \$sel_ref0]
}
set rmsf [measure rmsf \$sel_rms first $START_FRAME last -1]
for {set j 0} {\$j <  [\$sel_rms num] } {incr j} {
    puts \$outfile "[expr {\$j+1}] [lindex \$rmsf \$j]"
}
close \$outfile
EOF
((counter++))
done

cat >> tcl <<EOF
quit
EOF
vmd -dispdev text -e tcl

rm tcl
