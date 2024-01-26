
all: manual_5.3.pdf 

manual_5.3.pdf:manual_5.3.tex
	pdflatex  manual_5.3.tex
	bibtex manual_5.3
	pdflatex   manual_5.3.tex
	pdflatex  manual_5.3.tex
clean: 
	rm -f  *.dvi *.log *.aux *.toc *.bbl *.out *.blg