#!/bin/bash
#########################################
## Description: Main body of a "boss run" VMD/tcl measure script
## Author: 
##         Kevin C. Chan (work@skblnw.com) Feb 2022
## Usage: understand, modify and source it!
## Units: A
#########################################

# Selection for alpha-chain
SEL1="segname PROD and resid 94 to 102"
# Selection for beta-chain
SEL2="segname PROE and resid 94 to 102"
# Selection for MHC
SEL3="segname PROA"
# Selection for peptide
SEL4="segname PROC"
# Selection for TCR
SEL5="segname PROD PROE"

PDB="$1"
TRJ="$2"
OUTPUT="$3"
[ $# -ne 3 ] && { echo -e "mkvmd> Usage: $0 [PDB] [TRJ] [OUTPUT]\n       By default, the selections are:\n       Selection 1: $SEL1\n       Selection 2: $SEL2"; exit 1; }

if [ ! -f $PDB ]; then
    echo -e "$PDB \nStructure not found!"
    exit 0
fi

if [ ! -f $TRJ ]; then
    echo -e "$TRJ \nTrajectory not found!"
    exit 0
fi

rm ${OUTPUT}_a ${OUTPUT}_b ${OUTPUT}_mhc ${OUTPUT}_tcr
cat > tcl << EOF
proc countTCRPepContact { nn sel_input } {
    set sel [atomselect top "$SEL4"]
    puts "Selected [llength [lsort -unique -integer [\$sel get residue]]] residues"

    set nlist []
    set dlist {3.5 4 4.5}
    foreach ii [lsort -unique -integer [\$sel get residue]] {
        set tmplist 0
        foreach dd \$dlist {
            set selheavy [atomselect top "noh \$sel_input and same residue as within \$dd of residue \$ii" frame \$nn]
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

# Load packages for calculating principle axis (if needed)
package require Orient 
namespace import Orient::orient 

mol new $PDB waitfor all
mol addfile $TRJ waitfor all
set total_frame [molinfo top get numframes]

puts "mkvmd> Computing something..."
for {set nn 0} {\$nn < \$total_frame} {incr nn} {

    # /------------------------------------------------/
    # /     Where you really have to use your brain    /
    # /------------------------------------------------/

    set outf [open ${OUTPUT}_a "a"]
    set outf2 [open ${OUTPUT}_b "a"]
    set outf3 [open ${OUTPUT}_mhc "a"]
    set outf4 [open ${OUTPUT}_tcr "a"]
    # Write TIME at the very first of a line
    set out_line [format "%d" \$nn]
    # Definition of segments
    # This determines number of columns in one output file
    # Write to file
    puts \$outf "\$out_line [countTCRPepContact \$nn "$SEL1"]"
    puts \$outf2 "\$out_line [countTCRPepContact \$nn "$SEL2"]"
    puts \$outf3 "\$out_line [countTCRPepContact \$nn "$SEL3"]"
    puts \$outf4 "\$out_line [countTCRPepContact \$nn "$SEL5"]"
    # Remember to close the file
    close \$outf
    close \$outf2
    close \$outf3
    close \$outf4
}

quit
EOF

vmd -dispdev text -e tcl
rm tcl
