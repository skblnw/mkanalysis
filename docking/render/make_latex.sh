bash make_latex_figure.sh > fig
sed -n '/^\\documentclass/,/^%Add figures here$/p' template.tex > head
sed -n '/^%End figures here$/,/^\end{document}$/p' template.tex > tail
cat head fig tail > combine.tex
rm -f fig head tail

pdflatex combine.tex