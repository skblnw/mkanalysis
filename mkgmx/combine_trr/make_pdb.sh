GMX=/home/kevin/opt/gromacs-5.1.1-MPI-CUDA-single/bin/gmx_mpi

filename=../../../system.gro

$GMX editconf -f $filename -o initial.pdb -n ../../index.ndx
