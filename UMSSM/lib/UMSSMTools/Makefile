.PHONY: clean all

ifneq ($(MAKECMDGOALS),clean)
include  ../../../CalcHEP_src/FlagsForMake
endif

all: ../umhdecay

cLib = $(CALCHEP)/lib
../umhdecay:	umhdecay.F  sources/umhdecay.a ../aLib.a  $(cLib)/dynamic_me.a ../../work/work_aux.a  ../../../CalcHEP_src/lib/libSLHAplus.a ../../../sources/micromegas.a
	$(FC) $(FFLAGS) -o $@ $^ 

sources/umhdecay.a: sources/*.f sources/*.c
	$(MAKE) -C sources
	
clean:
	$(MAKE) -C sources clean
	rm -f  ../umhdecay
