.PHONY: clean all

ifneq ($(MAKECMDGOALS),clean)
include  ../../../CalcHEP_src/FlagsForMake
endif

all: ../umssmtoolslib

cLib = $(CALCHEP)/lib
../umssmtoolslib: umssmtools.F  sources/sources.a  $(cLib)/dynamic_me.a ../../work/work_aux.a ../aLib.a  ../../../CalcHEP_src/lib/libSLHAplus.a ../../../lib/micromegas.a ../../../CalcHEP_src/lib/ntools.a
	$(FC) $(FFLAGS) -o $@ $^ 

sources/sources.a: sources/*.f sources/*.c
	$(MAKE) -C sources
	
clean:
	$(MAKE) -C sources clean
	rm -f  ../umssmtoolslib
