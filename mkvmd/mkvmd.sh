#!/bin/bash
#########################################
## Description: Main body of a "boss run" VMD/tcl measure script
## Author: 
##         Kevin C. Chan (work@skblnw.com) Feb 2022
## Usage: understand, modify and source it!
## Units: A
#########################################

# Selection for 
SEL1=""
# Selection for 
SEL2=""
# Selection for
SEL3=""

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

rm ${OUTPUT}
cat > tcl << EOF
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
        lappend out_line []
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
