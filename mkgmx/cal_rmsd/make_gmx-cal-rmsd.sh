GMXPATH=

initial_gro=../../../system.gro
input_xtc=../combine_trr/md.xtc
index=../../index.ndx

gmx_mpi rms -s $initial_gro -f $input_xtc -n $index
