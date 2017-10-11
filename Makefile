
.PHONY: all clean flags

all:include/microPath.h
	$(MAKE) -C CalcHEP_src
	$(MAKE) -C sources

include/microPath.h:
	echo \#define micrO \"$(CURDIR)\"  > include/microPath.h
        
   
clean:  
	rm -f include/microPath.h
	./clean
	rm -rf */*.dSYM 

flags: 
	$(MAKE) -C CalcHEP_src flags