# Calc function as it is named
proc calc_diah {nn seg1 seg2 seg3 seg4 sel_input1 sel_input2 sel_input3 sel_input4} {
        set sel1 [atomselect top "segname A${seg1} and resid ${sel_input1} and name CA" frame $nn]
        set sel2 [atomselect top "segname A${seg2} and resid ${sel_input2} and name CA" frame $nn]
        set sel3 [atomselect top "segname A${seg3} and resid ${sel_input3} and name CA" frame $nn]
        set sel4 [atomselect top "segname A${seg4} and resid ${sel_input4} and name CA" frame $nn]
        
        set dih [dihedral [measure center $sel1] [measure center $sel2] [measure center $sel3] [measure center $sel4]]
        set res [format "%.2f" $dih]

        $sel1 delete
        $sel2 delete
        $sel3 delete
        $sel4 delete

        return $res
}