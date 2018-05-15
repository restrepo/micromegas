
all: manual_5.0.pdf 

manual_5.0.pdf:manual_5.0.tex
	latex  manual_5.0
	bibtex manual_5.0
	latex  manual_5.0
	latex  manual_5.0
	dvipdf  manual_5.0.dvi
clean: 
	rm -f  *.dvi *.log *.aux *.toc *.bbl *.out *.blg