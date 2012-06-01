.PHONY: clean all

Tools = NMSSMTools_2.3.1

ifneq ($(MAKECMDGOALS),clean)
include ../../CalcHEP_src/FlagsForMake
endif

aLibOBJ = $(notdir $(subst .c,.o,$(wildcard *.c))) fort.o


all:  aLib.a  nmhdecay nmspec 

aLib.a:aLib.a($(aLibOBJ))
	$(RANLIB) aLib.a
	rm -f $(aLibOBJ)

aLib.a(calchep_func.o): lpath.h

lpath.h:
	echo \#define LPATH \"`pwd`\" > lpath.h

nmhdecay: $(Tools)/nmhdecay.F $(Tools)/sources/*.f
	echo $(MAKE) -C  $(Tools) ../nmhdecay
	$(MAKE) -C  $(Tools) ../nmhdecay
 
nmspec: $(Tools)/nmspec.F $(Tools)/sources/*.f
	$(MAKE) -C  $(Tools) ../nmspec

clean: 
	rm -f lpath.h *.o *.a  *.exe  
	$(MAKE) -C $(Tools) clean
