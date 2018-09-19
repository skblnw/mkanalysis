count=0
echo -e '\\begin{figure*}\n'
echo -e '\\captionsetup[subfigure]{labelformat=empty}\n'
echo -e '\\setlength{\\lineskip}{1ex}\n'
cat namelist | while read line
do
    score=$(echo $line | awk '{print $2}')
    name=ZINC$line
    filename=ZINC$line.jpg

    let count+=1
    caption=`echo $name`
        echo -e '\\subfloat['$caption']{\\includegraphics[width=0.25\\textwidth]{'$filename'}}'
        echo -e '\\hspace{.125\\textwidth}'
    if (( $count % 18 == 0 )); then
        echo -e '\n\\end{figure*}\n'
        echo -e '\\begin{figure*}\n'
        echo -e '\\captionsetup[subfigure]{labelformat=empty}\n'
        echo -e '\\setlength{\\lineskip}{1ex}\n'
    fi
done
echo -e '\n\\end{figure*}\n'
