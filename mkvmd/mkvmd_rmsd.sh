#!/bin/bash
# Kev @ 2019
# Bash script to calculate RMSD profile

[ $# -eq 0 ] && { echo "mkvmd> Usage: $0 [PDB] [TRJ] [PREFIX] [REF]"; exit 1; }

PDB="$1"
TRJ="$2"
PREFIX="$3"
REF="$4"

# SELREF="noh segname PROA"
# REGIONS=("noh segname PROA" "noh segname PROB" "noh segname PROC" "noh segname PROA and resid 1 to 180")
SELREF="noh chain A"
REGIONS=("noh chain A" "noh chain B" "noh chain C" "noh chain A and resid 1 to 180")
FILENAMES=("proa" "prob" "proc" "proahelix")


[ ! -z "PREFIX" ] && echo "mkvmd> Please specify prefix for output files"

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
set sel_rmsd [atomselect 0 "${region}"]
set outfile [open $filename "w"]
puts \$outfile "# reference: $SELREF; calculated region: $region"
for {set i 0} {\$i<\$num_frames} {incr i} {
    \$sel_all frame \$i
    \$sel_ref frame \$i
    \$sel_rmsd frame \$i
    \$sel_all move [measure fit \$sel_ref \$sel_ref0]
    set rmsd [measure rmsd \$sel_rmsd \$sel_rmsd0 weight mass]
    puts \$outfile "\$i \t \$rmsd"
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
