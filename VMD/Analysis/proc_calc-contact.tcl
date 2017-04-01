proc calc_contact {nn seg1 seg2 sel_input1 sel_input2} {
    set sel1 [atomselect top "segname ${seg1} and resid ${sel_input1} and name CA" frame $nn]
    set sel2 [atomselect top "segname ${seg2} and resid ${sel_input2} and name CA" frame $nn]
    set rlist [measure contacts 10 $sel1 $sel2]
    
    set col1 [lindex $rlist 0]
    set col2 [lindex $rlist 1]

    set res [format "%d" [llength $col1]]
    return $res
}