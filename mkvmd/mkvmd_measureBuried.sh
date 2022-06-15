#!/bin/bash
#########################################
## Description: Main body of a "boss run" VMD/tcl measure script
## Author: 
##         Kevin C. Chan (work@skblnw.com) Feb 2022
## Usage: understand, modify and source it!
## Units: A
#########################################

SEL1="segname PEPT"
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

 # and not name CA C N O HN HA

rm ${OUTPUT}_area ${OUTPUT}_perc
cat > tcl << EOF
proc measureBuried { nn resid } {
    set rest [atomselect top "$SEL1 and resid \$resid and not name CA C N O HN HA HT1 HT2 HT3 OT1 OT2" frame \$nn]
    set selpep [atomselect top "$SEL1" frame \$nn]
    set selall [atomselect top "$SEL1 or $SEL2" frame \$nn]
    set percent_buried [format "%.2f" [expr 1 - [measure sasa 1.4 \$selall -restrict \$rest] / [measure sasa 1.4 \$selpep -restrict \$rest]]]
    set area_buried [format "%.2f" [expr [measure sasa 1.4 \$selpep -restrict \$rest] - [measure sasa 1.4 \$selall -restrict \$rest]]]
    # puts "ResID \$resid ([lsort -unique [\$rest get resname]]): Buried=\$area_buried"
    if { \$nn == 0 } {
        puts "ResID \$resid ([lsort -unique [\$rest get resname]]): Unbound=[format "%3.0f" [measure sasa 1.4 \$selpep -restrict \$rest]], Bound=[format "%3.0f" [measure sasa 1.4 \$selall -restrict \$rest]], Buried=[format "%3.0f" \$area_buried]"
    } else {
        puts -nonewline "[format "%3.0f" \$area_buried] + "
    }

    return [list \$percent_buried \$area_buried]
}

# /------------------/
# /     Main Body    /
# /------------------/

# !!!Important!!!
# Deleting existing files as we APPEND instead of trashing and opening new files
# Make sure you will delete all the existing files
# !!!Important!!!
# eval file delete [glob output/*.dat]
# set OUTPUT_DIR [exec date +%Y%m%d%H%M%S]
# exec mkdir -p $OUTPUT_DIR

# Load packages for calculating principle axis (if needed)
#package require Orient 
#namespace import Orient::orient 

mol new $PDB waitfor all
mol addfile $TRJ waitfor all
set total_frame [molinfo top get numframes]

puts "mkvmd> Computing something..."
for {set nn 0} {\$nn < \$total_frame} {incr nn} {
    set nframe [expr \$nn + 0]

    # /------------------------------------------------/
    # /     Where you really have to use your brain    /
    # /------------------------------------------------/

    # Reset index of selections
    set nsel 0
    
    # Definition of selections
    # You may use "...", e.g. "1 to 10", instead of one integer
    # This determines <number of output files>
    foreach {sel_input1} {1} {sel_input2} {1} {
      # Output file name
      # Increase index of selections by 1
      incr nsel
      set outf1 [open ${OUTPUT}_perc "a"]
      set outf2 [open ${OUTPUT}_area "a"]
      # Write TIME at the very first of a line
      set out_line1 [format "%d" \$nframe]
      set out_line2 [format "%d" \$nframe]
      # Definition of segments
      # You are free to use {seg1} {...} {seg2} {...}
    # This determines <number of columns in one output file>
      set total_area_buried 0
      foreach resid {$SEQUENCE} {
          # Call calc funtion you like
          set output [measureBuried \$nn \$resid]
          lappend out_line1 [lindex \$output 0]
          lappend out_line2 [lindex \$output 1]
          incr total_area_buried [format "%.0f" [lindex \$output 1]]
      }
      puts "Total Buried=\$total_area_buried"
      # Write to file
      puts \$outf1 "\$out_line1"
      puts \$outf2 "\$out_line2 \$total_area_buried"
      # Remember to close the file
      close \$outf1
      close \$outf2
    }
}

quit
EOF

vmd -dispdev text -e tcl
rm tcl
