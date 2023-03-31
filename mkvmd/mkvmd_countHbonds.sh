#!/bin/bash
#########################################
## Description: Main body of a "boss run" VMD/tcl measure script
## Author: 
##         Kev (contact@skblnw.com) Dec 2015
## Usage: understand, modify and bash it!
## Units: A
#########################################

SELTEXT1="segname PROC and resid 3 4 5 6 7"
SELTEXT2="name OW"

PDB="$1"
TRJ="$2"
OUTPUT="$3"
[ $# -ne 3 ] && { echo "mkvmd> Usage: $0 [PDB] [TRJ] [OUTPUT]"; echo "SEL1: $SELTEXT1; SEL2: $SELTEXT2"; exit 1; }

files=("$PDB" "$TRJ")
for file in "${files[@]}"; do
    if [ ! -f "$file" ]; then
        echo -e "$file \nStructure not found!"
        exit 1
    fi
done

rm $OUTPUT
cat > tcl << EOF
proc countHbonds {nn} {
    set sel1 [atomselect top "${SELTEXT1} and not backbone" frame \$nn]
    set sel2 [atomselect top "${SELTEXT2}" frame \$nn]
    set count1 [llength [lindex [measure hbonds 3.5 30 \$sel1 \$sel2] 0]]
    set count2 [llength [lindex [measure hbonds 3.5 30 \$sel2 \$sel1] 0]]
    \$sel1 delete
    \$sel2 delete

    return [expr \$count1 + \$count2]
}

# /------------------/
# /     Main Body    /
# /------------------/

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
rm tcl
