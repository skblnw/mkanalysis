mol new ../psfgen/d268g.psf waitfor all
mol addfile ../psfgen/d268g.pdb waitfor all
set sel1 [atomselect top all]

mol new ionized.psf waitfor all
mol addfile ionized.pdb waitfor all
set sel2 [atomselect top {segname "WT.*" ION}]

package require topotools
set mol [::TopoTools::selections2mol "$sel1 $sel2"]
set sel [atomselect $mol all]

animate write psf ionized-2.psf $mol
animate write pdb ionized-2.pdb $mol
