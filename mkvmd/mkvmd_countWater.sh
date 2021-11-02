#########################################
## Description: Main body of a "boss run" VMD/tcl measure script
## Author: 
##         Kev (contact@skblnw.com) Dec 2015
## Usage: understand, modify and bash it!
## Units: A
#########################################

PDB="$1"
TRJ="$2"
OUTPUT="$3"
[ $# -ne 3 ] && { echo "mkvmd> Usage: $0 [PDB] [TRJ] [OUTPUT]"; exit 1; }

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
proc countWater { nn sel_input } {
    set sel [atomselect top "\$sel_input"]
    puts "Selected [llength [lsort -unique -integer [\$sel get residue]]] residues"

    set nlist []
    set dlist {3 3.5 4}
    foreach ii [lsort -unique -integer [\$sel get residue]] {
        set tmplist 0
        foreach dd \$dlist {
            set selwater [atomselect top "name OH2 OW and within \$dd of residue \$ii" frame \$nn]
            incr tmplist [\$selwater num]
        }
        lappend nlist [expr \$tmplist / [llength \$dlist]]
    }
    
    set total 0
    foreach nxt \$nlist { incr total \$nxt }
    puts "Total number of water: \$total"
    return \$nlist
}

# /------------------/
# /     Main Body    /
# /------------------/

# /---------------------------------------------------------------------------------/
# /  Deleting existing files as we APPEND instead of trashing and opening new files /
# /---------------------------------------------------------------------------------/
# eval file delete [glob output/*.dat]
# set OUTPUT_DIR [exec date +%Y%m%d%H%M%S]
# exec mkdir -p \$OUTPUT_DIR

# Load packages for calculating principle axis if needed
# package require Orient 
# namespace import Orient::orient 

# Load your structure and frames
mol new $PDB waitfor all
mol addfile $TRJ waitfor all step 77
set total_frame [molinfo top get numframes]

for {set nn 0} {\$nn < \$total_frame} {incr nn} {
    set nframe [expr \$nn + 0]

    # /-------------------------------------------------/
    # /     Where you really have to use your brain     /
    # /-------------------------------------------------/
    # Uncomment these two lines if you need PA
    # set selref [atomselect top "protein and name CA" frame \$nn]
    # set Iref [Orient::calc_principalaxes \$selref]
    # /--------------------------------------------------------------/
    # /                       Atom Selections                        /
    # /   You may use "...", e.g. "1 to 10", instead of one integer  /
    # /          This determines <number of output columns>            /
    # /--------------------------------------------------------------/
    foreach {sel_input} {"segname MUT"} {
        set outf [open $OUTPUT "a"]
        # Write to file
        puts \$outf [countWater \$nn \$sel_input]
        # Remember to close the file
        close \$outf
    }
}
quit
EOF

vmd -dispdev text -e tcl
rm -f tcl
