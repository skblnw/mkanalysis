#!/bin/bash

DIR_WORK=$(pwd)
# e.g. /home/PHARMACY/chan.773/playground/cpf1/docking/vina/ligand/pbc_pdbqt/res
DIR_DOCKRES=
# e.g. /home/PHARMACY/chan.773/playground/cpf1/run_complex/mkanalysis/receptor/complex_280ns.pdb
RECEPTOR=
# e.g. _complex
SUFFIX=

cat > $DIR_WORK/cal${SUFFIX}.tcl << EOF

mol new combine${SUFFIX}.pdb
set reslist [lsort -unique [[atomselect top "same residue as {protein within 3 of resname UNL UNK}"] get resid]]
set ctype [[atomselect top "resid \$reslist and name CA"] get resname]
#set crna [lsort -unique [[atomselect top "same residue as {chain B and within 3 of resname UNL UNK}"] get resid]]
set of [open "contact${SUFFIX}" w]
puts \$of \$reslist
close \$of
set of [open "ctype${SUFFIX}" w]
puts \$of \$ctype
close \$of
set of [open "crna${SUFFIX}" w]
puts \$of \$crna
close \$of

EOF

mkdir -p combine
cd combine
for filename in $(cat $DIR_WORK/namelist)
do
    echo $filename
    mkdir -p $filename
    cd $filename
    sed -n '/^MODEL 1$/,/^ENDMDL$/p' $DIR_DOCKRES/$filename/out.pdbqt > tmp
    grep '^ATOM' tmp > tmp2
    mv tmp2 tmp

    head -n -1 $RECEPTOR > combine${SUFFIX}.pdb
    cat tmp >> combine${SUFFIX}.pdb
    echo -e "TER\nENDMDL" >> combine${SUFFIX}.pdb

    vmd -dispdev text -eofexit < $DIR_WORK/cal${SUFFIX}.tcl
    #mv combine${SUFFIX}.pdb combine.pdb
    rm -f tmp

    cd ..
done

