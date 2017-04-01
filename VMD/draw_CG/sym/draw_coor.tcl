set fp [open "com_1-11.dat" r]

set i 0
while {[gets $fp line] >= 0} {
	draw color red
	draw sphere $line radius 3
	set tmp [vecadd $line {4 4 4}]
	incr i
	draw text $tmp "P$i"
}
close $fp

set fp [open "com_12-20.dat" r]

while {[gets $fp line] >= 0} {
	draw color iceblue
	draw sphere $line radius 3
}
close $fp

set fp [open "com_21-31.dat" r]

while {[gets $fp line] >= 0} {
	draw color green
	draw sphere $line radius 3
}
close $fp

set fp [open "com_32-40.dat" r]

while {[gets $fp line] >= 0} {
	draw color cyan
	draw sphere $line radius 3
}
close $fp
