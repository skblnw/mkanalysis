GMXPATH=

initial_gro=../../../system.gro
input_xtc=../combine_trr/md.xtc
index=../../index.ndx
begin_ps=400

gmx_mpi rmsf -s $initial_gro -f $input_xtc -n $index -b $begin_ps -res
