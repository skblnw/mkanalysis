set psfname [lindex $argv 0]
set pdbname [lindex $argv 1]
set dcdname [lindex $argv 2]
set seltext [lindex $argv 3]
set outname [lindex $argv 4]

set molnum [mol new $psfname waitfor all]
mol addfile $pdbname waitfor all
mol addfile $dcdname waitfor all

#package require mdff
#mdff check -mol 0 -frames all -rmsd -rmsdseltext $seltext -refpdb $refname -rmsdfile ${outname}_rmsd.dat

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
