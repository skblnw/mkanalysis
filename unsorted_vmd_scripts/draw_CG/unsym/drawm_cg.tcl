# #########################################
# ## Description: TCL script to draw spheres according to some coordinates. Usually use for drawing CG sites.
# ## Author: Kevin May 2014
# ## Usage: source draw_coor.tcl
# ## Input: 3-column-matrix coordinate file, x y z
# ## Output: nice drawing of CG sites
# ## Units: 
# ## Other Notes: modify according to your needs
# #########################################


for {set start 1} {$start <= 1} {incr start} {
	draw delete all
	# draw materials on
	set n 1
	set i 1
	set c 6
	set x A
	set fp [open unsym+last.xyz r]
	set end [expr {$start+37}]
	set middle1 [expr {$start+1}]
	set middle2 [expr {$start+2}]
	while {[gets $fp line] >= 0} {
		# puts $line
		if {$n >=$start && $n<=$end} {
		  if {$c==2} {set c 10}
		  if {$c>=18  && $c%2==0} {incr c}
		  draw color $c
		  # incr c
		  set tmp [vecadd $line {-5 6 6}]
		  draw text $tmp P${n} size 2.5 thickness 5
		  draw sphere $line radius 2
		  
		  if {$i > 1} {draw cylinder $origin $line radius 0.3}
		  set origin $line
		  incr i
		  
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
		incr n
	}
	close $fp

# display resetview
# set string "P${start}+P${middle1}+P${middle2}+P${end}+dih-top.bmp"

# scale to 0.035

# render Tachyon tmp.dat [format {"C:\Program Files (x86)\University of Illinois\VMD\\tachyon_WIN32.exe"  -mediumshade -trans_max_surfaces 1 tmp.dat -res 2048 1024 -format BMP -o C:/Users/ukevi/Documents/Dropbox/Quick_Work_Directory/BAR-PH/draw_cg/40/%s} $string]

# rotate x to -90

# set string "P${start}+P${middle1}+P${middle2}+P${end}+dih-side.bmp"

# render Tachyon tmp.dat [format {"C:\Program Files (x86)\University of Illinois\VMD\\tachyon_WIN32.exe"  -mediumshade -trans_max_surfaces 1 tmp.dat -res 2048 1024 -format BMP -o C:/Users/ukevi/Documents/Dropbox/Quick_Work_Directory/BAR-PH/draw_cg/40/%s} $string]



# display resetview
# set string "P${start}+P${end}+bon-top.bmp"

# scale to 0.035

# render Tachyon tmp.dat [format {"C:\Program Files (x86)\University of Illinois\VMD\\tachyon_WIN32.exe"  -mediumshade -trans_max_surfaces 1 tmp.dat -res 2048 1024 -format BMP -o C:/Users/ukevi/Documents/Dropbox/Quick_Work_Directory/BAR-PH/draw_cg/40/%s} $string]

# rotate x to -90

# set string "P${start}+P${end}+bon-side.bmp"

# render Tachyon tmp.dat [format {"C:\Program Files (x86)\University of Illinois\VMD\\tachyon_WIN32.exe"  -mediumshade -trans_max_surfaces 1 tmp.dat -res 2048 1024 -format BMP -o C:/Users/ukevi/Documents/Dropbox/Quick_Work_Directory/BAR-PH/draw_cg/40/%s} $string]


}
