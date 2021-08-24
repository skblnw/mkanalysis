proc calEnergy { seg1 seg2 seltext1 seltext2 } {
    package require namdenergy
    set sel1 [atomselect top "segname P${seg1} and resid ${seltext1}"]
    set sel2 [atomselect top "segname P${seg2} and resid ${seltext2}"]

    namdenergy -elec -vdw -nonb -sel $sel1 $sel2 -ofile "output/ener-pair${ipair}-${ii}.dat" -switch 11 -cutoff 13 -par par_all36_prot.prm -exe /home/kevin/opt/NAMD_2.10_Linux-x86_64/namd2
}
