#!/bin/bash
#########################################
## Description: Main body of a "boss run" VMD/tcl measure script
## Author: 
##         Kev (contact@skblnw.com) Dec 2015
## Usage: understand, modify and bash it!
## Units: A
#########################################

SELTEXT1="segname PROA and resid 423 to 433"
SELTEXT2="segname PROP and resid 1286 to 1300"

PDB="$1"
TRJ="$2"
OUTPUT="$3"
[ $# -ne 3 ] && { echo "mkvmd> Usage: $0 [PDB] [TRJ] [OUTPUT]"; echo "SEL1: $SELTEXT1; SEL2: $SELTEXT2"; exit 1; }

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
proc lcompare { list1 list2 } {
    set outlist {}
    foreach ii \$list1 {
        if {[lsearch -exact \$list2 \$ii] != -1} {
            lappend outlist \$ii
        }
    }
    if {[llength \$outlist] != 0} {
        return \$outlist
    } else {
        return -1
    }
}

proc countHbonds {nn} {
    set outlist {}

    set sel1 [atomselect top "${SELTEXT1} and {backbone or name H}" frame \$nn]
    set sel2 [atomselect top "${SELTEXT2} and {backbone or name H}" frame \$nn]
    set count1 [llength [lindex [measure hbonds 3.5 30 \$sel1 \$sel2] 0]]
    set count2 [llength [lindex [measure hbonds 3.5 30 \$sel2 \$sel1] 0]]
    \$sel1 delete
    \$sel2 delete

    return [expr \$count1 + \$count2]
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
        puts \$outf "[format "%d" \$nframe] [countHbonds \$nn]"
        # Remember to close the file
        close \$outf
    }
}

quit
EOF

vmd -dispdev text -e tcl
rm -f tcl
