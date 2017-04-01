proc vmd_draw_arrow {mol start end} {
    # an arrow is made of a cylinder and a cone
    set middle [vecadd $start [vecscale 0.9 [vecsub $end $start]]]
    graphics $mol cylinder $start $middle radius 0.35
    graphics $mol cone $middle $end radius 0.55
}

set ff [open "pvec-1.dat" r]
set fc [open "com.dat" r]

while {[gets $ff vector]>=0 && [gets $fc coor]>=0} {
	set vector [vecscale $vector 120]
	set desn [vecadd $coor $vector]
	draw color white
	draw arrow $coor $desn
}