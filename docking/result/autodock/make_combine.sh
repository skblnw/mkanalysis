for filename in $(cat namelist)
do
    dir=$(dirname $filename)
    cd $dir
    sed -n '/^USER    Cluster Rank = 1$/,/^ENDMDL/p' 01.dlg > test
    grep '^ATOM' test > test2
    mv test2 test

    head -n -2 ../../../../complex/complex_280ns.pdb > combine.pdb
    cat test >> combine.pdb
    echo "TER\nENDMDL" >> combine.pdb

    rm -f test
    cp ../render.vmd .
    sh ../make_render_vmd.sh

    cd ..
done
