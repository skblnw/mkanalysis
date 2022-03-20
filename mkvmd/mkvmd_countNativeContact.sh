#!/bin/bash
#########################################
## Description: Main body of a "boss run" VMD/tcl measure script
## Author: 
##         Kev (contact@skblnw.com) Dec 2015
## Usage: understand, modify and bash it!
## Units: A
#########################################

SELTEXT1b="segname PROA and resid 436 to 449"
SELTEXT1="segname PROA and resid 423 to 434"
SELTEXT2="segname PROP and resid 1284 to 1300"
SELTEXT2b="segname PROP and resid 1046 to 1054"
SELTEXT2c="segname PROP and resid 1269 to 1284"

CUTOFF=5

PDB="$1"
TRJ="$2"
OUTPUT="$3"
REF="$4"
[ $# -ne 4 ] && { echo "mkvmd> Usage: $0 [PDB] [TRJ] [OUTPUT] [REF]"; exit 1; }

if [ ! -f $PDB ]; then
    echo -e "$PDB \nStructure not found!"
    exit 0
fi

if [ ! -f $TRJ ]; then
    echo -e "$TRJ \nTrajectory not found!"
    exit 0
fi

if [ ! -f $REF ]; then
    echo -e "$REF \nTrajectory not found!"
    exit 0
fi

rm $OUTPUT
cat > tcl << EOF

proc countNativeContact { nn sel1 sel2 ref } {
    set sel [atomselect top "noh \$sel2 and same residue as {within $CUTOFF of {\$sel1}}" frame \$nn]

    foreach ii [\$ref get resid] { set tarRes(\$ii) 1 }

    set count 0
    foreach res [\$sel get resid] {
        if [info exist tarRes(\$res)] {
            incr count
        }
    }
    
    set total [llength [ \$ref get resid ]]
    return [format "%4.2f" [expr {\$count / (\$total + 0.0)}]]
}

# /------------------/
# /     Main Body    /
# /------------------/

# /------------------------------------------------/
# / Deleting existing files as we APPEND instead of trashing and opening new files
# /------------------------------------------------/
# eval file delete [glob output/*.dat]
# set OUTPUT_DIR [exec date +%Y%m%d%H%M%S]
# exec mkdir -p \$OUTPUT_DIR

# Load packages for calculating principle axis if needed
# package require Orient 
# namespace import Orient::orient 

# Load your reference structure
mol new $REF waitfor all
set region1 [atomselect top "noh $SELTEXT1b and same residue as {within $CUTOFF of {$SELTEXT1}}"]
set region2a [atomselect top "noh $SELTEXT2 and same residue as {within $CUTOFF of {$SELTEXT1}}"]
set region2b [atomselect top "noh $SELTEXT2b and same residue as {within $CUTOFF of {$SELTEXT2}}"]
set region3 [atomselect top "noh $SELTEXT2c and same residue as {within $CUTOFF of {$SELTEXT2}}"]

# Load your structure and frames
mol new $PDB waitfor all
mol addfile $TRJ waitfor all
set total_frame [molinfo top get numframes]

puts "mkvmd> Computing something..."
for {set nn 0} {\$nn < \$total_frame} {incr nn} {
    set nframe [expr \$nn + 0]

    # /------------------------------------------------/
    # /     Where you really have to use your brain    /
    # /------------------------------------------------/
    # Uncomment these two lines if you need PA
    # set selref [atomselect top "protein and name CA" frame \$nn]
    # set Iref [Orient::calc_principalaxes \$selref]
    # /------------------------------------------------/
    # /                 Atom Selections
    # / You may use "...", e.g. "1 to 10", instead of one integer
    # / This determines <number of output files>
    # /------------------------------------------------/
    foreach {sel_input1} {1} {sel_input2} {1} {
        set outf [open $OUTPUT "a"]
        # Write to file
        puts \$outf "[format "%d" \$nframe] [countNativeContact \$nn "$SELTEXT1" "$SELTEXT1b" \$region1] [countNativeContact \$nn "$SELTEXT1" "$SELTEXT2" \$region2a] [countNativeContact \$nn "$SELTEXT2" "$SELTEXT2b" \$region2b] [countNativeContact \$nn "$SELTEXT2" "$SELTEXT2c" \$region3]"
        # Remember to close the file
        close \$outf
    }
}

quit
EOF

vmd -dispdev text -e tcl
rm -f tcl
