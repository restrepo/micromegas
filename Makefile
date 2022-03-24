
.PHONY: all clean flags


all:include/microPath.h
	$(MAKE) -C CalcHEP_src  MICROMEGAS=MICROMEGAS
	$(MAKE) -C sources

include/microPath.h:
	echo \#define micrO \"$(CURDIR)\"  > include/microPath.h
        
   
clean:  
	./clean
flags: 
	$(MAKE) -C CalcHEP_src flags