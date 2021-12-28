#########################################
## Description: Main body of a "boss run" VMD/tcl measure script
## Author: 
##         Kev (cchan2242-c@my.cityu.edu.hk) Dec 2015
## Usage: understand, modify and source it!
## Units: A
#########################################

SEL1="segname PROA"

PDB="$1"
TRJ="$2"
OUTPUT="$3"
[ $# -ne 3 ] && { echo -e "mkvmd> Usage: $0 [PDB] [TRJ] [OUTPUT]\n       By default, the selections are:\n       Selection 1: $SEL1"; exit 1; }

if [ ! -f $PDB ]; then
    echo -e "$PDB \nStructure not found!"
    exit 0
fi

if [ ! -f $TRJ ]; then
    echo -e "$TRJ \nTrajectory not found!"
    exit 0
fi

 # and not name CA C N O HN HA

rm ${OUTPUT}
cat > tcl << EOF
proc measureSasa { nn probe_size resid } {
    set rest [atomselect top "$SEL1 and resid \$resid" frame \$nn]
    set selprot [atomselect top "$SEL1" frame \$nn]
    set sasa [measure sasa \$probe_size \$selprot -restrict \$rest]
    if { \$sasa != 0 } {
        return 1
    } else {
        return 0
    }
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


set selprot [atomselect top "$SEL1"]
puts "mkvmd> These residues will be scanned: "
set tmplist []
foreach resid [lsort -unique -integer [\$selprot get resid]] { 
    set rest [atomselect top "$SEL1 and resid \$resid"]
    set sasa [measure sasa 5 \$selprot -restrict \$rest]
    if { \$sasa == 0 } { 
        lappend tmplist \$resid 
    }
}
set reslist []
foreach resid \$tmplist {  
    set rest [atomselect top "$SEL1 and resid \$resid"]
    set sasa [measure sasa 2.0 \$selprot -restrict \$rest]
    if { \$sasa != 0 } { 
        puts -nonewline "\$resid "
        lappend reslist \$resid 
    }
}
puts ""
puts "mkvmd> [llength \$reslist] residues in total"

puts "mkvmd> Computing something..."
for {set nn 0} {\$nn < \$total_frame} {incr nn} {
    set nframe [expr \$nn + 0]

    # /------------------------------------------------/
    # /     Where you really have to use your brain    /
    # /------------------------------------------------/
    
    # Definition of selections
    # You may use "...", e.g. "1 to 10", instead of one integer
    # This determines <number of output files>
    foreach {sel_input1} {1} {sel_input2} {1} {
      # Output file name
      set outf [open ${OUTPUT} "a"]
      # Write TIME at the very first of a line
      set out_line [format "%d" \$nframe]
      # Definition of segments
      # You are free to use {seg1} {...} {seg2} {...}
      # This determines <number of columns in one output file>

      if { \$nn == 0 } {
        foreach resid \$reslist {
            lappend out_line \$resid
        }
        puts \$outf "\$out_line"
        set out_line [format "%d" \$nframe]
      }

      foreach resid \$reslist {
        # Call calc funtion you like
        # set out_line [format "%.1f" \$probe_size]
        set count 0
        foreach probe_size {3.0 3.5 4.0} {
            incr count [measureSasa \$nn \$probe_size \$resid]
        }
        lappend out_line \$count
      }
      puts \$outf "\$out_line"

      # Write to file
      # Remember to close the file
      close \$outf
    }
}

quit
EOF

vmd -dispdev text -e tcl
rm tcl
