################################
## TCL script to get index of single chain protein
## Kevin Apr 2014
## Usage: vmd -dispdev text -e getind.tcl
## Output: list of index corresponding to atomselect options
################################

mol new .pdb type pdb waitfor all

set sel [atomselect top "segname and resid "]

set file [open index.txt w]
foreach indices [$sel get index] {
    puts $file "$indices"
}
close $file

#$sel writepdb initial.pdb

quit
