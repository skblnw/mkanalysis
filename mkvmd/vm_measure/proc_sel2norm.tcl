# Calc function as it is named
# 
# sel2 o   o sel3
#  vec1 \ / vec2
#        o
#       sel1
proc sel2norm {nn seg1 seg2 seg3 sel_input1 sel_input2 sel_input3} {
        set sel1 [atomselect top "segname P${seg1} and resid ${sel_input1} and name CA" frame $nn]
        set sel2 [atomselect top "segname P${seg2} and resid ${sel_input2} and name CA" frame $nn]
        set sel3 [atomselect top "segname P${seg3} and resid ${sel_input3} and name CA" frame $nn]
        
        set vec1 [vecsub [measure center $sel2 weight mass] [measure center $sel1 weight mass]]
        set vec2 [vecsub [measure center $sel3 weight mass] [measure center $sel1 weight mass]]
        set norm [veccross $vec1 $vec2]

        $sel1 delete
        $sel2 delete
        $sel3 delete

        return $norm
}
