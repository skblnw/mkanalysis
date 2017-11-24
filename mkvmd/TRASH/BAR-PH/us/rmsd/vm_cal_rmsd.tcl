set psfname [lindex $argv 0]
set dcdname [lindex $argv 1]
set seltext [lindex $argv 2]
set outname [lindex $argv 3]

set molnum [mol new $psfname waitfor all]
mol addfile $dcdname waitfor all

set num_frames [molinfo $molnum get numframes]
set a0 [atomselect top $seltext frame 0]
set b [atomselect top $seltext]
set outfile [open "${outname}_rmsd.dat" "w"]
for {set i 0} {$i<$num_frames} {incr i} {
    $b frame $i
    $b move [measure fit $b $a0]
    set rmsd1 [measure rmsd $b $a0]
    puts $outfile "$i \t $rmsd1 "
}
close $outfile

quit
