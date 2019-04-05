mol new ../initial/initial.pdb
mol addfile ../combine_dcd/cls1-gs3-water-protein-35ns-s50.dcd first 500 last -1 waitfor all
puts [molinfo top get numframes]
animate delete beg 0 end 0
puts [molinfo top get numframes]

set num_frames [molinfo top get numframes]
set a0 [atomselect top protein frame 0]
set b [atomselect top protein]
for {set ii 0} { $ii < $num_frames} {incr ii} {
    $b frame $ii
    $b move [measure fit $b $a0]
}

set sel [atomselect top "name CA"]
[atomselect top "name CA"] set {x y z} [measure avpos $sel]

$sel writepdb avg.pdb
