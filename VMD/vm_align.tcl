
set frames [molinfo top get numframes]
set all [atomselect top all]
set ref [atomselect top all frame 0]

# Align to first frame
if {1} {
puts "Info) Done wrapping, now align selection to the first frame"
for {set i 1} {$i < $frames} {incr i} {
  $all frame $i
  $all update

  $all move [measure fit $all $ref]
}
}

animate write trr out.trr

