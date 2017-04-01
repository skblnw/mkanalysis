for ii in 616 630 640 650 660 670 680 695 710 725 740 755 765 770 785 800 815 830 845 860 875 890 905 920 935 950 965 980
do
    vmd -dispdev text -e vm_make_ref_pdb.tcl -args ../../output/us-z${ii}-0.restart.coor ./output2/us-z$ii-0.pdb
done
