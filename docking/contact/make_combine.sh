#!/bin/bash

DIR_DOCKRES=/home/PHARMACY/chan.773/playground/cpf1/docking/vina/ligand/pbc_pdbqt/res_complex
RECEPTOR=/home/PHARMACY/chan.773/playground/cpf1/run_complex/mkanalysis/receptor/complex_280ns.pdb
DIR_RENDER=/home/PHARMACY/chan.773/playground/cpf1/docking/vina/ligand/pbc_pdbqt/render
DIR_OUTPUT=/home/PHARMACY/chan.773/Dropbox/QWD/CRISPR-Cpf1/cpf1/docking/ligand/pbc/render

cat > cal.tcl << EOF

mol new combine.pdb
set of [open "contact_complex" w]
puts \$of [lsort -unique [[atomselect top "same residue as {protein within 3 of resname UNL UNK}"] get resid]]
close \$of

EOF

for filename in $(cat namelist_complex)
do
    ZINC=$(basename $filename)
    echo $ZINC
    mkdir -p $ZINC
    cd $ZINC
    sed -n '/^MODEL 1$/,/^ENDMDL$/p' $DIR_DOCKRES/$ZINC/out.pdbqt > tmp
    grep '^ATOM' tmp > tmp2
    mv tmp2 tmp

    head -n -1 $RECEPTOR > combine.pdb
    cat tmp >> combine.pdb
    echo -e "TER\nENDMDL" >> combine.pdb

    vmd -dispdev text -eofexit < ../cal.tcl
    mv combine.pdb combine_complex.pdb
    rm -f tmp

    cd ..
done
