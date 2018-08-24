set sel [atomselect top "protein"]
set selcenter [atomselect top "protein"]
set centerX [measure center $selcenter]

for {set ii 0} {$ii < 360} {incr ii 30} {
	set matrixX [trans center $centerX axis x -10]
	$sel move $matrixX
}