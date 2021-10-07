#########################################
## Description: Main body of a "boss run" VMD/tcl measure script
## Author: 
##         Kev (cchan2242-c@my.cityu.edu.hk) Dec 2015
## Usage: understand, modify and source it!
## Units: A
#########################################

SEL1="segname PROA"
SEL2="segname MUT"
SEL_COMBINE="segname PROA MUT"

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

rm $OUTPUT
cat > tcl << EOF
proc measureSasa { nn } {
    set sel1 [atomselect top "$SEL1" frame \$nn]
    set sel2 [atomselect top "$SEL2" frame \$nn]
    set selall [atomselect top "$SEL_COMBINE" frame \$nn]
    
    return [format "%.2f" [expr -([measure sasa 1.4 \$selall] - [measure sasa 1.4 \$sel1] - [measure sasa 1.4 \$sel2]) / 2.0]]
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
      set outf [open $OUTPUT "a"]
      # Write TIME at the very first of a line
      set out_line [format "%d" \$nframe]
      # Definition of segments
      # You are free to use {seg1} {...} {seg2} {...}
    # This determines <number of columns in one output file>
      foreach {seg1 seg2} {1 1} {
          # Call calc funtion you like
          lappend out_line [measureSasa \$nn]
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
rm tcl
