mol new ../combine_dcd/initial.psf waitfor all
mol addfile ../combine_dcd/initial.pdb waitfor all
mol addfile ../combine_dcd/npt-pf100ps.dcd waitfor all

set sel [atomselect top protein]
$sel frame last
$sel writepdb final.pdb

quit