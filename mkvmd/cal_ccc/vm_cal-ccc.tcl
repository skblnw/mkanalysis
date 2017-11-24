set molnum [mol new ../combine_dcd/initial.psf waitfor all]
mol addfile ../combine_dcd/initial.pdb waitfor all
mol addfile ../combine_dcd/mdff-s50.dcd waitfor all

package require mdff

set num_frames [molinfo $molnum get numframes]
set sel [atomselect top "all"]
set outfile [open "ccc.dat" "w"]
for {set i 0} {$i<$num_frames} {incr i} {
    $sel frame $i
    set ccc1 [mdffi cc $sel -res 14 -i ../../../6tetramer_clip_rotated.mrc -thresholddensity 0.5]
    puts $outfile "$i \t $ccc1 "
}

quit
