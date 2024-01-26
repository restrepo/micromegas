
all: manual_5.3.pdf 

manual_5.3.pdf:manual_5.3.tex
	latex  manual_5.3
	bibtex manual_5.3
	latex  manual_5.3
	latex  manual_5.3
	dvipdf  manual_5.3.dvi
clean: 
	rm -f  *.dvi *.log *.aux *.toc *.bbl *.out *.blg