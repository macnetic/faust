all: Faust.pdf

pdf: Faust.pdf

Faust.pdf: clean Faust.tex
	pdflatex Faust.tex
	makeindex Faust.idx
	pdflatex Faust.tex
	latex_count=8 ; \
	while egrep -s 'Rerun (LaTeX|to get cross-references right)' Faust.log && [ $$latex_count -gt 0 ] ;\
	    do \
	      echo "Rerunning latex...." ;\
	      pdflatex Faust.tex ;\
	      latex_count=`expr $$latex_count - 1` ;\
	    done
	makeindex Faust.idx
	pdflatex Faust.tex


clean:
	rm -f *.ps *.dvi *.aux *.toc *.idx *.ind *.ilg *.log *.out *.brf *.blg *.bbl Faust.pdf
