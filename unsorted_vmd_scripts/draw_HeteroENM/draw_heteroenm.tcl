#draw delete all

set MAP [open "map.xyz" r]
while {[gets $MAP line] > 0} {
	set n [lindex $line 0]
	set comx($n) [lindex $line 1]
	set comy($n) [lindex $line 2]
	set comz($n) [lindex $line 3]
}
close $MAP


set IN [open "cgk.dat" r]
while {[gets $IN line] >= 0} {
	set i [lindex $line 0]
	set j [lindex $line 1]
	set bl [lindex $line 2]
	set k [lindex $line 3]
	
	if {$k > 0.1} {
	set start [list $comx($i) $comy($i) $comz($i)]
	set end [list $comx($j) $comy($j) $comz($j)]
	set mid [vecadd $start $end]
	set mid [vecscale 0.5 $mid]
	if {$k > 2.0} {draw color red
	} elseif {$k > 1.0} {draw color orange
	} elseif {$k > 0.1} {draw color lime
	} elseif {$k > 0.01} {draw color iceblue
	} else {draw color white}
	draw cylinder $start $end radius 0.2
	#draw text $mid "$k" size 0.8 thickness 1.5
    }

}
close $IN