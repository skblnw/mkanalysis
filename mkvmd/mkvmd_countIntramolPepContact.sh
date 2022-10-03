#!/bin/bash
#########################################
## Description: Main body of a "boss run" VMD/tcl measure script
## Author: 
##         Kevin C. Chan (work@skblnw.com) Feb 2022
## Usage: understand, modify and source it!
## Units: A
#########################################

SEL1="segname PROC"
SEL2="segname PROA PROB"

PDB="$1"
TRJ="$2"
OUTPUT="$3"
[ $# -ne 3 ] && { echo -e "mkvmd> Usage: $0 [PDB] [TRJ] [OUTPUT]\n       By default, the selections are:\n       Selection 1: $SEL1\n       Selection 2: $SEL2\n       Selection combined: $SEL_COMBINE"; exit 1; }

if [ ! -f $PDB ]; then
    echo -e "$PDB \nStructure not found!"
    exit 0
fi

if [ ! -f $TRJ ]; then
    echo -e "$TRJ \nTrajectory not found!"
    exit 0
fi

 # and not name CA C N O HN HA

rm ${OUTPUT}_heavy ${OUTPUT}_snp ${OUTPUT}_sp
cat > tcl << EOF
proc countIntramolContact { nn sel_input } {
    set sel [atomselect top "$SEL1"]
    puts "Selected [llength [lsort -unique -integer [\$sel get residue]]] residues"

    set nlist []
    set dlist {4}
    foreach ii [lsort -unique -integer [\$sel get residue]] {
        set tmplist 0
        foreach dd \$dlist {
            set selheavy [atomselect top "\$sel_input and not residue \$ii and within \$dd of {residue \$ii and sidechain}" frame \$nn]
            incr tmplist [\$selheavy num]
        }
        lappend nlist [expr \$tmplist / [llength \$dlist]]
    }
    
    set total 0
    foreach nxt \$nlist { incr total \$nxt }
    puts "Total number of nearby residues: \$total"
    return \$nlist
}

# /------------------/
# /     Main Body    /
# /------------------/

mol new $PDB waitfor all
mol addfile $TRJ waitfor all
set total_frame [molinfo top get numframes]

puts "mkvmd> Computing something..."
for {set nn 0} {\$nn < \$total_frame} {incr nn} {

    # /------------------------------------------------/
    # /     Where you really have to use your brain    /
    # /------------------------------------------------/

    set outf [open ${OUTPUT}_heavy "a"]
    set outf2 [open ${OUTPUT}_snp "a"]
    set outf3 [open ${OUTPUT}_sp "a"]
    # Write TIME at the very first of a line
    set out_line [format "%d" \$nn]
    # Definition of segments
    # This determines number of columns in one output file
    # Write to file
    puts \$outf "\$out_line [countIntramolContact \$nn "noh $SEL1 and sidechain"]"
    puts \$outf2 "\$out_line [countIntramolContact \$nn "$SEL1 and sidechain and name \"C.*\""]"
    puts \$outf3 "\$out_line [countIntramolContact \$nn "$SEL1 and sidechain and name \"O.*\" \"N.*\""]"
    # Remember to close the file
    close \$outf
    close \$outf2
    close \$outf3
}

quit
EOF

vmd -dispdev text -e tcl
rm tcl
