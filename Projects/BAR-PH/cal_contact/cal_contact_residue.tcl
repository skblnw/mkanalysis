set molnum [mol new ../initial/initial.pdb waitfor all]
mol addfile ../combine_dcd/actin-npt1-193-s50.dcd waitfor all

set nframe [molinfo $molnum get numframes]

for {set i 0} {$i < $nframe} {incr i} {
    puts "Frame $i"
    foreach {seg} {3 4 5 6 7 8 9 10 11 12 13} {
        set sel1 [atomselect top "segname A${seg} and resid 40 and name CA"]
        set sel2 [atomselect top "not segname A${seg} and name CA"]
        $sel1 frame $i
        $sel2 frame $i
        
        set res [measure contacts 10 $sel1 $sel2]
        
        set col1 [lindex $res 0]
        set col2 [lindex $res 1]
        
        set outf [open ./output/contact-A$seg-$i.dat "w"]
        puts "Writing results to ./output/contact-A$seg-$i.dat"
        for {set j 0} {$j < [llength $col1]} {incr j} {
            set idx1 [lindex $col1 $j]
            set idx2 [lindex $col2 $j]
            
            set sel3 [atomselect top "index $idx1"]
            set out1 [$sel3 get resid]
            set out3 [$sel3 get segname]
            set sel4 [atomselect top "index $idx2"]
            set out2 [$sel4 get resid]
            set out4 [$sel4 get segname]
            
            puts $outf "$out1\t$out2\t$out3\t$out4"
        }
        close $outf
    }
}
quit
