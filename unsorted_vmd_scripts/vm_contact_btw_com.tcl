set pos1 [measure center $sel1]
set pos2 [measure center $sel2]
set dist [veclength [vecsub $pos2 $pos1]]
if {$dist < $cutoff} { incr contacts; } 
