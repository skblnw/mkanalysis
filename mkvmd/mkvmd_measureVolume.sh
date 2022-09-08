#########################################
## Description: Measure volume of an HLA groove by solvating and subtracting water molecules from the complex. 
## Author: 
##         Kev (cchan@zju.edu.cn) Sep 2022
## Usage: understand, modify and source it!
## Units: A
#########################################

PEPT="segname PROC"

PDB="$1"
TRJ="$2"
OUTPUT="$3"
[ $# -ne 3 ] && { echo -e "mkvmd> Usage: $0 [PDB] [TRJ] [OUTPUT]\n       By default, the peptide selections are:\n       $SEL"; exit 1; }

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
proc measureVolume { nn } {
    set sel [atomselect 0 protein frame \$nn]
    set mol [::TopoTools::selections2mol "\$sel"]
    animate write psf bound.psf \$mol
    animate write pdb bound.pdb \$mol
    set sel1 [atomselect 0 "protein and not $PEPT" frame \$nn]
    set sel2 [atomselect 0 "protein and $PEPT and not within 5 of {segname PROA PROB}" frame \$nn]
    set mol [::TopoTools::selections2mol "\$sel1 \$sel2"]
    animate write psf free.psf \$mol
    animate write pdb free.pdb \$mol

    set dlist {2.3 2.4 2.5}
    set vol 0
    foreach dd \$dlist {
        solvate free.psf free.pdb -t 1 -o solvated_free -rotate -b \$dd
        solvate bound.psf bound.pdb -t 1 -o solvated_bound -rotate -b \$dd
        mol new solvated_free.pdb
        set sel1 [atomselect top "name OH2"]
        mol new solvated_bound.pdb
        set sel2 [atomselect top "name OH2"]
        incr vol [expr ([\$sel1 num] - [\$sel2 num]) * 30]
    }
    
    return [expr \$vol / [llength \$dlist]]
}

# /------------------/
# /     Main Body    /
# /------------------/

package require topotools
package require solvate

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
          lappend out_line [measureVolume \$nn]
      }
      # Write to file
      puts "measureVolume> \$out_line"
      puts \$outf "\$out_line"
      # Remember to close the file
      close \$outf
    }
}

quit
EOF

vmd -dispdev text -e tcl | grep "^measureVolume>"
rm tcl
