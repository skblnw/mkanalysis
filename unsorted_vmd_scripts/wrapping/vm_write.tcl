set seltext {protein or segname "L.*" "W." MEMB ION or {segname "WT.*" and name OH2 and within 10 of protein}}
set outname "noW"

mol new ../../ionized.psf
mol addfile ../../ionized.pdb

set sel [atomselect top $seltext]

package require topotools
set mol1 [::TopoTools::selections2mol "$sel"]
animate write psf initial-$outname.psf $mol1
animate write pdb initial-$outname.pdb $mol1

set file [open index-$outname.txt w]
foreach indices [$sel get index] {
  puts $file "$indices"
}

close $file
quit
