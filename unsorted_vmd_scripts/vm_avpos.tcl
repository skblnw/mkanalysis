mol new .psf
mol addfile .dcd first 500 last -1 waitfor all
puts [molinfo top get numframes]
#animate delete beg 0 end 0
puts [molinfo top get numframes]

set sel0 [atomselect top protein frame 0]
set sel [atomselect top protein]
$sel0 set {x y z} [measure avpos $sel]

$sel0 writepdb avg.pdb