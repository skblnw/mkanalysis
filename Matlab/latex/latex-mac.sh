count=0
echo  '\\begin{figure}[htbp]\n'
for filename in *.png
do
    let count+=1
    caption=`echo $filename | cut -d. -f1 | cut -d- -f2,3`
#    echo  '\\begin{figure}[htbp]\n\includegraphics[width=.3\\textwidth]{PLOT_subunit-front-sel'$ii'.png}\n\\caption{}\n\label{fig:}\n\\end{figure}\n'
        echo  '\\begin{subfigure}[b]{.3\\textwidth}\n\includegraphics[width=\\textwidth]{'$filename'}\n\\caption{'$caption'}\n\label{fig:}\n\\end{subfigure}'
    if (( $count % 12 == 0 )); then
        echo  '\n\\end{figure}\n'
        echo  '\\begin{figure}[htbp]\n'
    fi
done
echo  '\n\\end{figure}\n'