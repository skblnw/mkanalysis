################################
## TCL script to get index of single chain protein
## Kevin Apr 2014
## Usage: vmd -dispdev text -e getind.tcl
## Output: list of index corresponding to atomselect options
################################

mol new ../first/initial.pdb type pdb waitfor all

set sel [atomselect top "segname P1 P2 P7 P8 and resid 331 to 338 or segname P5 P6 P3 P4 and resid 230 to 249"]

set file [open output/index-I2-sel6.txt w]
foreach indices [$sel get index] {
    puts $file "$indices"
}
close $file

$sel writepdb output/initial-I2-sel6.pdb

quit
