
.PHONY: all clean flags

all:
	$(MAKE) -C CalcHEP_src
	$(MAKE) -C sources
   
clean:  
	rm -f manual24.log manual24.aux manual24.dvi manual24.ps manual24.pdf
	./clean

flags: 
	$(MAKE) -C CalcHEP_src flags