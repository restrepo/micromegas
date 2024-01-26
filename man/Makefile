
all: manual_6.0.pdf 

manual_6.0.pdf:manual_6.0.tex
	pdflatex  manual_6.0.tex
	bibtex manual_6.0
	pdflatex   manual_6.0.tex
	pdflatex   manual_6.0.tex
clean: 
	rm -f  *.dvi *.log *.aux *.toc *.bbl *.out *.blg