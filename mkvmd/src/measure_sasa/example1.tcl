mol new uc63/t22.pdb
foreach ii {25} {mol addfile uc63/t$ii.pdb}
foreach ii {10 13 14 15 16} {mol addfile uc64/t$ii.pdb}

set seln [atomselect top "protein and resid 1 to 197"]
set selc [atomselect top "protein and resid 198 to 300"]
set selp [atomselect top "protein"]

measure sasa 1.4 $seln -points pp
draw color red
foreach p $pp {draw point $p}
measure sasa 1.4 $selp -restrict $seln -points pp
draw color blue
foreach p $pp {draw point $p}