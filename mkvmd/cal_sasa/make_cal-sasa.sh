# -args <psf> <dcd> <outname>
vmd -dispdev text -e vm_cal-sasa.tcl -args ../combine_dcd/initial.pdb ../combine_dcd/npt-pf100ps.dcd sasa_npt
