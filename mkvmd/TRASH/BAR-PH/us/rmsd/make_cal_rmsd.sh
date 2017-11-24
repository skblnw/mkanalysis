module load vmd
for ii in 616 630 640 650 660 670 680 695 710 725 740 755 765 770 785 800 815 830 845 860 875 890 905 920 935 950 965 980
do
vmd -dispdev text -e vm_cal_rmsd.tcl -args ../combine_dcd/output2/us-z$ii-0.pdb ../combine_dcd/output/us-z$ii-s100.dcd protein ./output/us-z${ii}
done
