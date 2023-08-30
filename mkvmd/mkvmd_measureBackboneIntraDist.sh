#!/bin/bash
#########################################
## Description: Main body of a "boss run" VMD/tcl measure script
## Author: 
##         Kevin C. Chan (work@skblnw.com) Feb 2022
## Usage: understand, modify and source it!
## Units: A
#########################################

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

echo "" > ${OUTPUT}.csv
cat > tcl << EOF

proc measureBackboneIntraDist { nn backbone_residue } {
    set distance_list {}
    set next_backbone_residue {}

    # These are resid of nitrogens
    foreach element \$backbone_residue {lappend next_backbone_residue [expr {\$element + 4}]} 

    set id1 [[atomselect top "resid \$backbone_residue and name O" frame \$nn] get index]
    set id2 [[atomselect top "resid \$next_backbone_residue and name N" frame \$nn] get index]
    foreach aa \$id1 bb \$id2 {
        set distance [measure bond [list \$aa \$bb] frame \$nn]
        lappend distance_list \$distance
    }
    return \$distance_list
}

# /------------------/
# /     Main Body    /
# /------------------/

# Load packages for calculating principle axis (if needed)
package require Orient 
namespace import Orient::orient 

# Load trajectory
mol new $PDB waitfor all
mol addfile $TRJ waitfor all
set total_frame [molinfo top get numframes]

puts "mkvmd> Computing something..."
for {set nn 0} {\$nn < \$total_frame} {incr nn} {

    # Define output filename
    set outf [open ${OUTPUT}.csv "a"]

    # These are resid of oxygens
    # set backbone_residue {310 311 312 313 314 315 316 317 318 319 320} 
    set backbone_residue {351 352 353 354 355 356 357 358 359 360 361} 
    set output [measureBackboneIntraDist \$nn \$backbone_residue]
    
    foreach element \$output {puts -nonewline \$outf [format "%.2f," \$element]}
    puts \$outf ""

    # Remember to close the file
    close \$outf

}

quit
EOF

vmd -dispdev text -e tcl
rm tcl
