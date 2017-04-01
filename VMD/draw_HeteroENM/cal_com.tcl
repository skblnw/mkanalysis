proc center_of_mass {selection} {
        # some error checking
        if {[$selection num] <= 0} {
                error "center_of_mass: needs a selection with atoms"
        }
        # set the center of mass to 0
        set com [veczero]
        # set the total mass to 0
        set mass 0
        # [$selection get {x y z}] returns the coordinates {x y z} 
        # [$selection get {mass}] returns the masses
        # so the following says "for each pair of {coordinates} and masses,
	#  do the computation ..."
        foreach coord [$selection get {x y z}] m [$selection get mass] {
           # sum of the masses
           set mass [expr $mass + $m]
           # sum up the product of mass and coordinate
           set com [vecadd $com [vecscale $m $coord]]
        }
        # and scale by the inverse of the number of atoms
        if {$mass == 0} {
                error "center_of_mass: total mass is zero"
        }
        # The "1.0" can't be "1", since otherwise integer division is done
        return [vecscale [expr 1.0/$mass] $com]
}

set bdindexs(1) 1
set bdindexs(2) 101
set bdindexs(3) 435
set bdindexs(4) 832
set bdindexs(5) 1165
set bdindexs(6) 1481
set bdindexs(7) 1886
set bdindexs(8) 2230
set bdindexs(9) 2540
set bdindexs(10) 2949
set bdindexs(11) 3353
set bdindexs(12) 3619
set bdindexs(13) 3931
set bdindexs(14) 4170
set bdindexs(15) 4274
set bdindexs(16) 4509
set bdindexs(17) 4767
set bdindexs(18) 4922
set bdindexs(19) 5091
set bdindexs(20) 5264
set bdindexs(21) 5461
set bdindexs(22) 5582
set bdindexs(23) 5904
set bdindexs(24) 6264
set bdindexs(25) 6610
set bdindexs(26) 6926
set bdindexs(27) 7305
set bdindexs(28) 7742
set bdindexs(29) 8083
set bdindexs(30) 8403
set bdindexs(31) 8728
set bdindexs(32) 9107
set bdindexs(33) 9392
set bdindexs(34) 9685
set bdindexs(35) 10033
set bdindexs(36) 10144
set bdindexs(37) 10372
set bdindexs(38) 10616
set bdindexs(39) 10807
set bdindexs(40) 10954
set bdindexs(41) 11142
set bdindexs(42) 11313
set bdindexs(43) 11445

set bdindexe(1) 100
set bdindexe(2) 434
set bdindexe(3) 831
set bdindexe(4) 1164
set bdindexe(5) 1480
set bdindexe(6) 1885
set bdindexe(7) 2229
set bdindexe(8) 2539
set bdindexe(9) 2948
set bdindexe(10) 3352
set bdindexe(11) 3618
set bdindexe(12) 3930
set bdindexe(13) 4169
set bdindexe(14) 4273
set bdindexe(15) 4508
set bdindexe(16) 4766
set bdindexe(17) 4921
set bdindexe(18) 5090
set bdindexe(19) 5263
set bdindexe(20) 5460
set bdindexe(21) 5581
set bdindexe(22) 5903
set bdindexe(23) 6263
set bdindexe(24) 6609
set bdindexe(25) 6925
set bdindexe(26) 7304
set bdindexe(27) 7741
set bdindexe(28) 8082
set bdindexe(29) 8402
set bdindexe(30) 8727
set bdindexe(31) 9106
set bdindexe(32) 9391
set bdindexe(33) 9684
set bdindexe(34) 10032
set bdindexe(35) 10143
set bdindexe(36) 10371
set bdindexe(37) 10615
set bdindexe(38) 10806
set bdindexe(39) 10953
set bdindexe(40) 11141
set bdindexe(41) 11312
set bdindexe(42) 11444
set bdindexe(43) 11725

set fo [open "map.xyz" w]
for {set j 1} {$j<=43} {incr j} {
    set sel [atomselect top "serial $bdindexs($j) to $bdindexe($j)"]
    $sel set beta [expr $j*100]
    puts $fo "$j [center_of_mass $sel]"
}
close $fo

set fp [open "map.xyz" r]