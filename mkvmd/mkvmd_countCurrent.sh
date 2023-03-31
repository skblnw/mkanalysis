#!/bin/bash
#########################################
## Description: Main body of a "boss run" VMD/tcl measure script
## Author: 
##         Kev (contact@skblnw.com) Dec 2015
## Usage: understand, modify and bash it!
## Units: A
#########################################

SELTEXT1b="segname PROA and resid 436 to 449"

CUTOFF=5

PDB="$1"
TRJ="$2"
OUTPUT="$3"
[ $# -ne 3 ] && { echo "mkvmd> Usage: $0 [PDB] [TRJ] [OUTPUT]"; exit 1; }

if [ ! -f $PDB ]; then
    echo -e "$PDB \nStructure not found!"
    exit 0
fi

if [ ! -f $TRJ ]; then
    echo -e "$TRJ \nTrajectory not found!"
    exit 0
fi

rm $OUTPUT
cat > tcl << EOF

proc countCurrency { nn jj } {
    set region_old [atomselect top "noh $SELTEXT1b and same residue as {within $CUTOFF of {$SELTEXT1}}" frame \$jj]
    set region_new [atomselect top "noh $SELTEXT1b and same residue as {within $CUTOFF of {$SELTEXT1}}" frame \$nn]

    foreach ii [\$region_old get resid] { set tarRes(\$ii) 1 }
    set count 0
    foreach res [\$region_new get resid] {
        if { ![info exist tarRes(\$res)] } {
            puts \$res
            puts [\$res get {x y z}]
        }
    }

    set total [llength [ \$ref get resid ]]
    return [format "%4.2f" [expr {\$count / (\$total + 0.0)}]]
}

# /------------------/
# /     Main Body    /
# /------------------/

# Load your structure and frames
mol new $PDB waitfor all
mol addfile $TRJ waitfor all
set total_frame [molinfo top get numframes]

puts "mkvmd> Computing something..."
for {set nn 1} {\$nn < \$total_frame} {incr nn} {
    set nframe [expr \$nn + 0]

    # /------------------------------------------------/
    # /     Where you really have to use your brain    /
    # /------------------------------------------------/
    # /------------------------------------------------/
    # /                 Atom Selections
    # / You may use "...", e.g. "1 to 10", instead of one integer
    # / This determines <number of output files>
    # /------------------------------------------------/

    set jj [expr \$nn - 1]

    set outf [open $OUTPUT "a"]
    # Write to file
    puts \$outf "[format "%d" \$nframe] [countCurrency \$nn \$jj]"
    # Remember to close the file
    close \$outf
}

quit
EOF

vmd -dispdev text -e tcl
rm -f tcl
