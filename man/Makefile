
all: manual_4.3.4.pdf 

manual_4.3.4.pdf:manual_4.3.4.tex
	latex  manual_4.3.4
	bibtex manual_4.3.4
	latex  manual_4.3.4
	latex  manual_4.3.4
	dvipdf  manual_4.3.4.dvi
clean: 
	rm -f  *.dvi *.log *.aux *.toc *.bbl *.out *.blg