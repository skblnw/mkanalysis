cat > cal.tcl << EOF

mol new combine.pdb
set of [open "out.contact" w]
puts \$of [lsort -unique [[atomselect top "same residue as {protein within 3 of chain d}"] get resid]]
close \$of

EOF

for filename in $(cat namelist)
do

    dir=$(dirname $filename)

    cd $dir
    vmd -dispdev text -eofexit < ../cal.tcl
    cd ..

done