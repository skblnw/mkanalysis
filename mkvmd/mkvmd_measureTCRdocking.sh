#!/bin/bash
#########################################
## Description: Main body of a "boss run" VMD/tcl measure script
## Author: 
##         Kev (cchan2242-c@my.cityu.edu.hk) Dec 2015
## Usage: understand, modify and source it!
## Units: A
#########################################

SEL1="segname ANTI"
SEL2="segname PROA PROB"
SEQUENCE="1 2 3 4 5 6 7 8 9"

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

rm ${OUTPUT}
cat > tcl << EOF
proc angle { a b } {
  # get Pi 
  global M_PI

   # Angle between two vectors 
  set cosine [expr [vecdot \$a \$b] / ( [veclength \$a] * [veclength \$b])]
  return [expr acos(\$cosine)*(180.0/\$M_PI)]
}

proc measureTCRDocking {nn} {
  set sel1 [atomselect top "segname PROD and resid 24 90" frame \$nn]
  set sel2 [atomselect top "segname PROE and resid 25 93" frame \$nn]
  set tcr [vecsub [measure center \$sel1 weight mass] [measure center \$sel2 weight mass]]

  set sel3 [atomselect top "segname PROA and resid 50 to 86 140 to 176" frame \$nn]
  set mhc [lindex [Orient::calc_principalaxes \$sel3] 2]

  \$sel1 delete
  \$sel2 delete
  \$sel3 delete
  return [format "%.2f" [angle \$tcr \$mhc]]
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
    set nframe [expr \$nn + 0]

    # /------------------------------------------------/
    # /     Where you really have to use your brain    /
    # /------------------------------------------------/

    set outf [open ${OUTPUT} "a"]
    # Write TIME at the very first of a line
    set out_line [format "%d" \$nframe]
    # Definition of segments
    # This determines number of columns in one output file
    set total_area_buried 0
    foreach {seg1 seg2} {1 1} {
        # Call calc funtion you like
        lappend out_line [measureTCRDocking \$nn]
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
