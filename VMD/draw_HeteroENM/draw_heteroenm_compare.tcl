#draw delete all

set MAP [open "map.xyz" r]
while {[gets $MAP line] > 0} {
	set n [lindex $line 0]
	set comx($n) [lindex $line 1]
	set comy($n) [lindex $line 2]
	set comz($n) [lindex $line 3]
}
close $MAP


set IN [open "water_sym.dat" r]
while {[gets $IN line] >= 0} {
	set i [lindex $line 0]
	set j [lindex $line 1]
	set ks [lindex $line 4]
	set ku [lindex $line 5]
	
	set start [list $comx($i) $comy($i) $comz($i)]
	set end [list $comx($j) $comy($j) $comz($j)]
	set mid [vecadd $start $end]
	set mid [vecscale 0.5 $mid]
	if {$ks > $ku} {draw color blue
	} else {draw color red}
	draw cylinder $start $end radius 0.3
	draw text $mid "$ks||$ku" size 0.8 thickness 1.5

}
close $IN