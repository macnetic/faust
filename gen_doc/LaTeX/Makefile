all: getting-started-with-faust-toolbox.pdf

pdf: getting-started-with-faust-toolbox.pdf

getting-started-with-faust-toolbox.pdf: clean getting-started-with-faust-toolbox.tex
	pdflatex getting-started-with-faust-toolbox.tex
	#makeindex getting-started-with-faust-toolbox.idx
	pdflatex getting-started-with-faust-toolbox.tex
	latex_count=8 ; \
	while egrep -s 'Rerun (LaTeX|to get cross-references right)' getting-started-with-faust-toolbox.log && [ $$latex_count -gt 0 ] ;\
	    do \
	      echo "Rerunning latex...." ;\
	      pdflatex getting-started-with-faust-toolbox.tex ;\
	      latex_count=`expr $$latex_count - 1` ;\
	    done
	#makeindex Faust.idx
	pdflatex getting-started-with-faust-toolbox.tex


clean:
	rm -f *.ps *.dvi *.aux *.toc *.idx *.ind *.ilg *.log *.out *.brf *.blg *.bbl getting-started-with-faust-toolbox.pdf
