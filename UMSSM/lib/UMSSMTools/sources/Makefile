.PHONY: clean all

ifneq ($(MAKECMDGOALS),clean)
include ../../../../CalcHEP_src/FlagsForMake
endif

all: sources.a

SRCF=$(wildcard *.f)
OBJF=$(SRCF:.f=.o)
SRCC=$(wildcard *.c)
OBJC=$(SRCC:.c=.o)
OBJ=$(OBJC) $(OBJF)

sources.a:sources.a($(OBJ))
	$(RANLIB) sources.a
	rm -f *.o
	
	
clean:
	rm -f *.a  *.o
