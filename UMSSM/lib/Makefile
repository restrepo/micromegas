.PHONY: clean all

Tools = UMSSMTools

all: aLib.a umssmtoolslib

ifneq ($(MAKECMDGOALS),clean)
include  ../../CalcHEP_src/FlagsForMake
endif

aLibOBJ = $(notdir $(subst .c,.o,$(wildcard *.c)))


aLib.a:aLib.a($(aLibOBJ))
	$(RANLIB) aLib.a
	rm -f $(aLibOBJ)

aLib.a(mo_func.o): lpath.h

lpath.h:
	echo \#define LPATH \"`pwd`\" > lpath.h

umssmtoolslib: $(Tools)/umssmtools.F $(Tools)/sources/*.f $(Tools)/sources/*.c
	echo $(MAKE) -C  $(Tools) ../umssmtoolslib
	$(MAKE) -C  $(Tools) ../umssmtoolslib

clean:
	rm -f lpath.h *.o *.a  *.exe
	rm -rf *.dSYM
	$(MAKE) -C $(Tools) clean
