set selp [atomselect top "protein"]
set selcenter [atomselect top "protein"]
set centerX [list [$selcenter get x] [$selcenter get y] [$selcenter get z]]

for {set ii 0} {x < 360} {incr ii 30} {
	set matrixX [trans center $centerX axis x -10]
	$sel move $matrixX
}