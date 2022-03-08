#!/bin/bash
#########################################
## Description: Main body of a "boss run" VMD/tcl measure script
## Author: 
##         Kev (cchan2242-c@my.cityu.edu.hk) Dec 2015
## Usage: understand, modify and source it!
## Units: A
#########################################

# Selections for disulfite in alpha (SEL1) and beta (SEL2) chain of TCR
SEL1="segname PROD and resid 89 to 93 and resname CYS and name CA"
SEL2="segname PROE and resid 89 to 93 and resname CYS and name CA"
# Selection for MHC alpha-helices
# Conventonally residues A:50 to 86 140 to 176 for MHC-I
# Conventonally residues A46 to 78, B:54 to 64 67 to 91 for MHC-II
# SEL4 is the peptide
SEL3="segname PROA and resid 50 to 86 140 to 176"
# SEL3="segname PROA and resid 46 to 78 or segname PROB and resid 54 to 64 67 to 91"
SEL4="segname PROC"

PDB="$1"
TRJ="$2"
OUTPUT="$3"
[ $# -ne 3 ] && { echo -e "mkvmd> Usage: $0 [PDB] [TRJ] [OUTPUT]\n       By default, the selections are:\n       Selection 1: $SEL1\n       Selection 2: $SEL2\n       Selection 3: $SEL3"; exit 1; }

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
  set sel1 [atomselect top "$SEL1" frame \$nn]
  set sel2 [atomselect top "$SEL2" frame \$nn]
  set tcr [vecsub [measure center \$sel2 weight mass] [measure center \$sel1 weight mass]]

  set sel3 [atomselect top "$SEL3" frame \$nn]
  set mhc [lindex [Orient::calc_principalaxes \$sel3] 2]

  set sel4 [atomselect top "$SEL4" frame \$nn]
  set vec [vecsub [lindex [\$sel4 get {x y z}] end] [lindex [\$sel4 get {x y z}] 0]]
  if { [angle \$mhc \$vec] > 90 } {
    set mhc [vecinvert \$mhc]
  }

  \$sel1 delete
  \$sel2 delete
  \$sel3 delete
  \$sel4 delete
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
