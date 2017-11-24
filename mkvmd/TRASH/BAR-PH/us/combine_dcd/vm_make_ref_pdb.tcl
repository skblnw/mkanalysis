set refname [lindex $argv 0]
set outname [lindex $argv 1]

mol new ../../ionized.psf waitfor all
mol addfile $refname waitfor all
set sel [atomselect top "protein"]
$sel writepdb $outname.pdb
quit
