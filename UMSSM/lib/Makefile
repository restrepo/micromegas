.PHONY: clean all

Tools = UMSSMTools

all: aLib.a umhdecay

ifneq ($(MAKECMDGOALS),clean)
include  ../../CalcHEP_src/FlagsForMake
endif

aLibOBJ = $(notdir $(subst .c,.o,$(wildcard *.c)))


aLib.a:aLib.a($(aLibOBJ))
	$(RANLIB) aLib.a
	rm -f $(aLibOBJ)

aLib.a(calchep_func.o): lpath.h

lpath.h:
	echo \#define LPATH \"`pwd`\" > lpath.h

umhdecay: $(Tools)/umhdecay.F $(Tools)/sources/*.f $(Tools)/sources/*.c
	echo $(MAKE) -C  $(Tools) ../umhdecay
	$(MAKE) -C  $(Tools) ../umhdecay

clean:
	rm -f lpath.h *.o *.a  *.exe  
	$(MAKE) -C $(Tools) clean
