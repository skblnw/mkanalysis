set psfname [lindex $argv 0]
set dcdname [lindex $argv 1]
set outname [lindex $argv 2]

set molnum [mol new $psfname waitfor all]
mol addfile $dcdname waitfor all

set num_frames [molinfo top get numframes]
set prot [atomselect top protein]
set surf [atomselect top "segname P1 and resid 137 to 170 27 to 35 60 to 75 and sidechain"]
set outfile [open "${outname}.dat" "w"]
for {set i 0} {$i<$num_frames} {incr i} {
    $prot frame $i
    $prot update
    $surf frame $i
    $surf update
    set sasa [measure sasa 1.40 $prot -restrict $surf -samples 100]
    puts $outfile "$i \t $sasa "
}
close $outfile

quit
