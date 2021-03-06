set sel1 [atomselect 0 "protein"]
set sel2 [atomselect 2 "not {water or ion}"]

package require topotools

set mol [::TopoTools::selections2mol "$sel1 $sel2"]
set selall [atomselect $mol all]

animate write psf system.psf $mol
animate write pdb system.pdb $mol

mol delete all
mol new system.psf waitfor all
mol addfile system.pdb waitfor all