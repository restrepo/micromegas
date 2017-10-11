  
.PHONY: all clean

ifneq ($(MAKECMDGOALS),clean)
include ../CalcHEP_src/FlagsForMake
endif

CC+= -I../include
all:../lib/micromegas.a 

OBJ = $(notdir $(subst .c,.o,$(filter-out $(special),$(wildcard *.c))))


#ifneq ($(FC),)
#   override OBJ += faux.o
#endif
../lib/micromegas.a(cs22.o):data/et_tab.inc

../lib/micromegas.a:../lib/micromegas.a($(OBJ)) 
	$(RANLIB) ../lib/micromegas.a 
	rm -f *.o

clean:
	rm -f ../lib/micromegas.a  *.o 
