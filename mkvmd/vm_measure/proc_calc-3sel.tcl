# Calc function as it is named
proc calc_angle {nn seg1 seg2 seg3 sel_input1 sel_input2 sel_input3} {
        set sel1 [atomselect top "segname A${seg1} and resid ${sel_input1} and name CA" frame $nn]
        set sel2 [atomselect top "segname A${seg2} and resid ${sel_input2} and name CA" frame $nn]
        set sel3 [atomselect top "segname A${seg3} and resid ${sel_input3} and name CA" frame $nn]
        
        set vec1 [vecsub [measure center $sel1] [measure center $sel2]]
        set vec2 [vecsub [measure center $sel3] [measure center $sel2]]
        set ang [angle $vec1 $vec2]
        set res [format "%.2f" $ang]

        $sel1 delete
        $sel2 delete
        $sel3 delete

        return $res
}
