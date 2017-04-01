set sel1 [atomselect 0 all] 
set sel2 [atomselect 1 all]

package require topotools 

set mol [::TopoTools::selections2mol "$sel1 $sel2"]
set sel [atomselect $mol all] 

animate write psf merged.psf $mol
animate write pdb merged.pdb $mol