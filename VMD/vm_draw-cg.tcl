# #########################################
# ## Description: TCL script to draw spheres according to some coordinates. Usually use for drawing CG sites.
# ## Author: Kevin May 2014
# ## Usage: source draw_coor.tcl
# ## Input: 3-column-matrix coordinate file, x y z
# ## Output: nice drawing of CG sites
# ## Units: 
# ## Other Notes: modify according to your needs
# #########################################

draw delete all
for {set start 1} {$start <= 1} {incr start} {
	draw delete all
	draw materials on
	draw material AOChalky
	set ii 0
	set c 0
	set cmap {"red" "blue" "cyan" "magenta"}
	set x A
	set fp [open com.dat r]
	set end [expr {$start+37}]
	set middle1 [expr {$start+1}]
	set middle2 [expr {$start+2}]
	while {[gets $fp line] >= 0} {
		# puts $line
		if {$n >=$start && $n<=$end} {
		  if {$c==2} {set c 10}
		  if {$c>=18  && $c%2==0} {incr c}
		  draw color [lindex $cmap $ii]
		  incr c
		  set tmp [vecadd $line {10 2 6}]
		  # draw text $tmp R${n} size 2.5 thickness 8
		  draw sphere $line radius 2
		  
		  draw color black
		  if {$ii > 0} {draw cylinder $origin $line radius 0.3}
		  set origin $line
		  
		  # if {$n==$start} {draw material GlassBubble;draw sphere $line radius 30}
		  } else {
			  draw color 6
			  draw sphere $line radius 2
			  # incr i
			  # draw text $tmp PB${i} size 2.5 thickness 1.5
			}

		# draw sphere $line radius 2
		
		# if {$n == 21} {
		  # set i 1
		  # set c 0
		  # set x B
		# }
		incr ii
	}
	close $fp
}
