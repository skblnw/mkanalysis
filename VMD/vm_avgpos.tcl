mol new ../initial/initial.pdb
mol addfile ../combine_dcd/cls1-gs3-water-protein-35ns-s50.dcd first 500 last -1 waitfor all
puts [molinfo top get numframes]
animate delete beg 0 end 0
puts [molinfo top get numframes]

set sel [atomselect top protein]
[atomselect top protein] set {x y z} [measure avpos $sel]

$sel writepdb avg.pdb
