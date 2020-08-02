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

set seltext(2) "resid 5 to 33 80 to 147 334 to 349"
set seltext(1) "resid 34 to 39 52 to 59"
set seltext(3) "resid 148 to 179 273 to 333"
set seltext(4) "resid 180 to 219 252 to 262"

set fo [open "com.dat" w]
for {set j 1} {$j<=4} {incr j} {
    set sel [atomselect top "$seltext($j)"]
    puts $fo [center_of_mass $sel]
}
close $fo