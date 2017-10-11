.PHONY: clean all

ifneq ($(MAKECMDGOALS),clean)
include ../../../../CalcHEP_src/FlagsForMake
endif

all: umhdecay.a

SRCF=$(wildcard *.f)
OBJF=$(SRCF:.f=.o)
SRCC=$(wildcard *.c)
OBJC=$(SRCC:.c=.o)
OBJ=$(OBJC) $(OBJF)

umhdecay.a:umhdecay.a($(OBJ))
	$(RANLIB) umhdecay.a
	rm -f *.o
	
	
clean:
	rm -f *.a  *.o
