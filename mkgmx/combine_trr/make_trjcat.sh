GMX=/home/kevin/opt/gromacs-5.1.1-MPI-CUDA-single/bin/gmx_mpi

TRJ_DIR=/share/data/kevin/martini/systems/Ups/dimer/run/output/

#filename="$TRJ_DIR/eq.xtc"

for ii in {1..5}
do
    filename="$filename $TRJ_DIR/md-$ii.xtc"
done

$GMX trjcat -f $filename -o md.xtc -n ../../index.ndx -dt 1000 -settime
