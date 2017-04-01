set molnum [mol new ../combine_dcd/cls2-4dimer.psf waitfor all]
mol addfile ../combine_dcd/cls2-gs3-water-protein-8ns-s50.dcd waitfor all

package require mdff

#mdff check -mol 0 -frames all -ccc -map ../combine_dcd/cls2.mrc -res 14 -cccfile cls2-4dimer-vac-gs3-result_ccc_global.dat
mdff check -mol 0 -frames all -ccc -map ../combine_dcd/cls2_updated.mrc -res 14 -cccfile cls2-4dimer-vac-gs3-result_ccc_local.dat -threshold 0.5

quit
