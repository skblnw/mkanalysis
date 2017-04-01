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
set bdindexs(2) 401
set bdindexs(3) 491
set bdindexs(4) 880
set bdindexs(5) 1615
set bdindexs(6) 2109
set bdindexs(7) 2373
set bdindexs(8) 2559
set bdindexs(9) 2724
set bdindexs(10) 3046
set bdindexs(11) 3467
set bdindexs(12) 3907
set bdindexs(13) 4274
set bdindexs(14) 4509
set bdindexs(15) 4767
set bdindexs(16) 4922
set bdindexs(17) 5091
set bdindexs(18) 5264
set bdindexs(19) 5461
set bdindexs(20) 5582
set bdindexs(21) 5904
set bdindexs(22) 6314
set bdindexs(23) 6396
set bdindexs(24) 6793
set bdindexs(25) 7519
set bdindexs(26) 8013
set bdindexs(27) 8294
set bdindexs(28) 8470
set bdindexs(29) 8628
set bdindexs(30) 8961
set bdindexs(31) 9392
set bdindexs(32) 9827
set bdindexs(33) 10144
set bdindexs(34) 10372
set bdindexs(35) 10616
set bdindexs(36) 10807
set bdindexs(37) 10954
set bdindexs(38) 11142
set bdindexs(39) 11313
set bdindexs(40) 11445

set bdindexe(1) 400
set bdindexe(2) 490
set bdindexe(3) 879
set bdindexe(4) 1614
set bdindexe(5) 2108
set bdindexe(6) 2372
set bdindexe(7) 2558
set bdindexe(8) 2723
set bdindexe(9) 3045
set bdindexe(10) 3466
set bdindexe(11) 3906
set bdindexe(12) 4273
set bdindexe(13) 4508
set bdindexe(14) 4766
set bdindexe(15) 4921
set bdindexe(16) 5090
set bdindexe(17) 5263
set bdindexe(18) 5460
set bdindexe(19) 5581
set bdindexe(20) 5903
set bdindexe(21) 6313
set bdindexe(22) 6395
set bdindexe(23) 6792
set bdindexe(24) 7518
set bdindexe(25) 8012
set bdindexe(26) 8293
set bdindexe(27) 8469
set bdindexe(28) 8627
set bdindexe(29) 8960
set bdindexe(30) 9391
set bdindexe(31) 9826
set bdindexe(32) 10143
set bdindexe(33) 10371
set bdindexe(34) 10615
set bdindexe(35) 10806
set bdindexe(36) 10953
set bdindexe(37) 11141
set bdindexe(38) 11312
set bdindexe(39) 11444
set bdindexe(40) 11725

set fo [open "com_1-11.dat" w]
for {set j 1} {$j<=11} {incr j} {
    set sel [atomselect top "serial $bdindexs($j) to $bdindexe($j)"]
    puts $fo [center_of_mass $sel]
}
close $fo

set fo [open "com_12-20.dat" w]
for {set j 12} {$j<=20} {incr j} {
    set sel [atomselect top "serial $bdindexs($j) to $bdindexe($j)"]
    puts $fo [center_of_mass $sel]
}
close $fo

set fo [open "com_21-31.dat" w]
for {set j 21} {$j<=31} {incr j} {
    set sel [atomselect top "serial $bdindexs($j) to $bdindexe($j)"]
    puts $fo [center_of_mass $sel]
}
close $fo

set fo [open "com_32-40.dat" w]
for {set j 32} {$j<=40} {incr j} {
    set sel [atomselect top "serial $bdindexs($j) to $bdindexe($j)"]
    puts $fo [center_of_mass $sel]
}
close $fo