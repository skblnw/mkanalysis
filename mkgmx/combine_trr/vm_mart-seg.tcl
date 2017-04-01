set sellist [list "1 to 410" "411 to 603"]

mol new initial.pdb waitfor all

set ii 0
foreach seltext $sellist {
	incr ii
	set sel [atomselect top "serial $seltext"]
	$sel set segname [join "P$ii"]
    $sel delete
}

set sel [atomselect top all]
$sel writepdb initial-seg.pdb