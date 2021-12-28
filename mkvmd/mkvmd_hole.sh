#!/bin/bash

[ $# -ne 3 ] && { echo "mkhole> Usage: $0 <pdb> <dcd> <output directory>"; exit 1; }
CATDCD='/home/kevin/opt/vmd/plugins/LINUXAMD64/bin/catdcd5.2/catdcd'
HOLE='hole'


cat > radius <<'EOF'
remark: Input file for program Tooshort.
remark: Contains bond and vdw radius records for
remark: for pdb atoms.
remark: format:
remark:   BOND C??? 0.8 
remark:   = bond radius of atom with first character C
remark:   is 0.8 angs.  ? is a wildcard character
remark: keyword VDWR for van der Waals radii
remark:    VDWR CA   GLY 1.925
remark: N.B. PUT SPECIFIC RECORDS BEFORE GENERAL
remark:
remark: van der Waals radii: AMBER united atom
remark: from Weiner et al. (1984), JACS, vol 106 pp765-768
BOND C??? 0.85
BOND N??? 0.75
BOND O??? 0.7
BOND S??? 1.1
BOND H??? 0.5
BOND P??? 1.0
remark: van der Waals radii
remark: all cb's are type C2 except thr, val, ile
VDWR CB   ALA 2.00
VDWR CB   THR 1.85
VDWR CB   VAL 1.85
VDWR CB   ILE 1.85
VDWR CB   ??? 1.925
remark: other C2 atoms (carbon with two aliphatic hydrogens)
VDWR CA   GLY 1.925
VDWR CG   GLU 1.925
VDWR CG   LYS 1.925
VDWR CD   LYS 1.925
VDWR CE   LYS 1.925
VDWR CG   PRO 1.925
VDWR CD   PRO 1.925
VDWR CG   MET 1.925
VDWR CG1  ILE 1.925
VDWR CG   GLN 1.925
VDWR CG   ARG 1.925
VDWR CD   ARG 1.925
remark: C3 atoms (carbon with two aliphatic hydrogens)
VDWR CD?  LEU 2.00
VDWR CE   MET 2.00
VDWR CG2  THR 2.00
VDWR CD1  ILE 2.00
VDWR CG2  ILE 2.00
VDWR CG?  VAL 2.00
remark: NH3 atom type N3
VDWR NZ   LYS 1.85
remark: sp2 oxygen atom type O
VDWR OE1  GLN 1.60
VDWR O    ??? 1.60
remark: acid oxygens O2
VDWR OE?  GLU 1.60
VDWR OD?  ASP 1.60
remark: general last
VDWR C??? ??? 1.85
VDWR O??? ??? 1.65
VDWR S??? ??? 2.00
VDWR N??? ??? 1.75
VDWR H??? ??? 1.00
VDWR P??? ??? 2.10
EOF

cat > tcl <<EOF
mol new $1 type pdb waitfor all
set sel [atomselect top "protein"]
set file [open index w]
foreach indices [\$sel get index] {puts \$file "\$indices"}
close \$file
quit
EOF
vmd -dispdev text -e tcl > /dev/null
rm tcl

cat > tcl <<'EOF'
mol new tmp.pdb
set sel [atomselect top all]
$sel writepdb tmp.pdb
quit
EOF

cat > hole.inp << EOF
coord tmp.pdb
radius radius
cvect 0.0 0.0 1.0
sphpdb out.sph
endrad 10
EOF

$CATDCD $2 > tmp
nframe=`grep "Total frames:" tmp | awk '{print $3}'`

[ -d $3 ] && { mv $3 BAK.$3; mkdir $3; } || mkdir $3
for ii in $(seq 1 $nframe)
do
    rm -f tmp.pdb
    echo -ne "\r> Doing HOLE computation for frame ${ii}/$nframe"
    $CATDCD -otype pdb -o tmp.pdb -i index -first $ii -last $ii -s $1 $2 > /dev/null
    vmd -dispdev text -e tcl > /dev/null 2>&1 
    $HOLE <hole.inp> out.txt
    egrep "mid-|sampled" out.txt | awk '{print $1" "$2}' > formatted.dat
    mv formatted.dat $3/$ii.dat
done

rm tmp.pdb out.sph out.sph.old out.txt 
echo "> Done."
