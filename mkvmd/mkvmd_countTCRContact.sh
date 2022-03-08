#!/bin/bash
#########################################
## Description: Main body of a "boss run" VMD/tcl measure script
## Author: 
##         Kevin C. Chan (work@skblnw.com) Feb 2022
## Usage: understand, modify and source it!
## Units: A
#########################################

# Selection for both alpha- and beta-chain
SEL1="segname PROD PROE"
# Selection for alpha-chain
SEL2="segname PROD"
# Selection for beta-chain
SEL3="segname PROE"
# Selection for MHC+peptide
SEL="segname PROA PROB PROC"

PDB="$1"
TRJ="$2"
OUTPUT="$3"
[ $# -ne 3 ] && { echo -e "mkvmd> Usage: $0 [PDB] [TRJ] [OUTPUT]\n       By default, the selections are:\n       TCR Selection: $SEL1\n       MHC Selection: $SEL"; exit 1; }

if [ ! -f $PDB ]; then
    echo -e "$PDB \nStructure not found!"
    exit 0
fi

if [ ! -f $TRJ ]; then
    echo -e "$TRJ \nTrajectory not found!"
    exit 0
fi

rm ${OUTPUT}
cat > tcl << EOF
proc countTCRContact { nn sel_input } {
    set dlist {3.5 4 4.5}
    set count 0
    foreach dd \$dlist {
        set selheavy [atomselect top "noh and $SEL and same residue as within \$dd of \$sel_input" frame \$nn]
        incr count [\$selheavy num]
        \$selheavy delete
    }
    lappend avgcount [expr \$count / [llength \$dlist]]
    
    return \$avgcount
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

    set outf [open ${OUTPUT} "a"]
    # Write TIME at the very first of a line
    set out_line [format "%d" \$nn]
    # Definition of segments
    # This determines number of columns in one output file
    foreach {seg1 seg2} {1 1} {
        # Call calc funtion you like
        lappend out_line [countTCRContact \$nn "$SEL1"] [countTCRContact \$nn "$SEL2"] [countTCRContact \$nn "$SEL3"]
    }
    # Write to file
    puts \$outf "\$out_line"
    # Remember to close the file
    close \$outf
}

quit
EOF

vmd -dispdev text -e tcl
rm tcl
