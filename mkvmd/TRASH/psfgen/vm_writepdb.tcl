#################################################
## TCL script to write a pdb of initial structure
## Kevin Apr 2014
## Usage: vmd -dispdev text -e <tcl>
## Output: a pdb file
#################################################

mol new ../../../ionized.pdb
set sel [atomselect top protein]
$sel writepdb initial.pdb
mol delete all
