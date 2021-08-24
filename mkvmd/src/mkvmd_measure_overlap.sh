#########################################
## Description: Main body of a "boss run" VMD/tcl measure script
## Author: 
##         Kev (cchan2242-c@my.cityu.edu.hk) Dec 2015
## Usage: understand, modify and source it!
## Units: A
#########################################

PDB="$1"
TRJ="$2"
REF="$3"
OUTPUT="$4"
[ $# -ne 4 ] && { echo "mkvmd> Usage: $0 [PDB] [TRJ] [REF] [OUTPUT]"; exit 1; }

if [ ! -f $PDB ]; then
    echo -e "$PDB \nStructure not found!"
    exit 0
fi

if [ ! -f $TRJ ]; then
    echo -e "$TRJ \nTrajectory not found!"
    exit 0
fi

if [ ! -f $REF ]; then
    echo -e "$REF \nStructure not found!"
    exit 0
fi


cat > tcl << EOF
source proc_findOverlap.tcl

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

# Load your structure and frames
mol new $REF
set ref [atomselect top "name CA and same residue as {chain B and within 6 of chain A}"]
set tt [llength [\$ref get resid]]

mol new $PDB waitfor all
mol addfile $TRJ waitfor all
set total_frame [molinfo top get numframes]

for {set nn 0} {\$nn < \$total_frame} {incr nn} {
    set nframe [expr \$nn + 0]

    # /------------------------------------------------/
    # /     Where you really have to use your brain    /
    # /------------------------------------------------/
    # Uncomment if you need PA
#    set selref [atomselect top "protein and name CA" frame \$nn]
#    set Iref [Orient::calc_principalaxes \$selref]

    # Reset index of selections
    set nsel 0
    
    # Definition of selections
    # You may use "...", e.g. "1 to 10", instead of one integer
    # This determines <number of output files>
    foreach {sel_input1} {1} {sel_input2} {1} {
      # Output file name
      # Increase index of selections by 1
      incr nsel
      set outf [open $OUTPUT "a"]
      # Write TIME at the very first of a line
      set out_line [format "%d" \$nframe]
      # Definition of segments
      # You are free to use {seg1} {...} {seg2} {...}
    # This determines <number of columns in one output file>
      foreach {seg1 seg2} {1 1} {
          # Call calc funtion you like
          lappend out_line [findOverlap \$nn \$ref 6 \$tt]
      }
      # Write to file
      puts \$outf "\$out_line"
      # Remember to close the file
      close \$outf
    }
}

quit
EOF

vmd -dispdev text -e tcl
rm -f tcl
