for filename in $(cat namelist)
do
    dir=$(dirname $filename)
    cp $dir/out.jpg jpg/$dir.jpg
done
