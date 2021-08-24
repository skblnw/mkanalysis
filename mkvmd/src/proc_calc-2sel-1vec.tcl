source proc_angle.tcl

# Calc function of angles between a vector (defined by sel1, sel2) and the principle axis (input Ip3)
proc calc_vec_pa {nn seg1 seg2 sel_input1 sel_input2 Ip3} {
        set sel1 [atomselect top "segname A${seg1} and resid ${sel_input1} and name CA" frame $nn]
        set sel2 [atomselect top "segname A${seg2} and resid ${sel_input2} and name CA" frame $nn]

        set vec [vecsub [measure center $sel1] [measure center $sel2]]
        set ang [angle $vec $Ip3]
        set res [format "%.2f" $ang]

        $sel1 delete
        $sel2 delete

        return $res
}
