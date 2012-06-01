.PHONY: clean all

ifneq ($(MAKECMDGOALS),clean)
include ../../../../CalcHEP_src/FlagsForMake
endif

all:	nmhdecay.a 

SRC = $(wildcard *.f)
OBJ = $(SRC:.f=.o)

nmhdecay.a:nmhdecay.a($(OBJ)) 
	$(RANLIB) nmhdecay.a 
	rm -f *.o

clean:
	rm -f *.a  *.o
