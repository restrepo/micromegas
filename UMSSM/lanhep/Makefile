#	This Makefile and the UMSSM model files which it compiles are used in LanHEP version 3.2.0
#	File varver.c of LanHEP version 3.2.0 has to be replaced by the one given in this directory to generate the same number of abbreviations in func1.mdl (thanks to Andrei Semenov).
#	If you want to compile it and use the resulting .mdl files in the UMSSM directory of micromegas, don't forget to :
#	1) Add '0' in column 'aux', for the 3 RH-neutrinos in file prtcls1.mdl;
#	2) Add the line '%Local!|' after the line defining the function 'RtoN6C2' in func1.mdl, to send all necessary functions to UMSSMTools.

FORTLIB = -lgfortran
SHLIBS = extfunc.so

.PHONY: all clean distclean

all: $(SHLIBS) ../../Packages/LanHEP/lhep
	../../Packages/LanHEP/lhep startup.src -ca -evl 2
#	add option "-allvrt" to also get 2-legs vertices

platform=$(shell uname)
ifeq ($(platform),Darwin)
extfunc.so: extfunc.c
	gcc -shared -o extfunc.so extfunc.c
else ifeq ($(platform),Linux)
extfunc.so: extfunc.c
	gcc -fPIC -shared -o extfunc.so extfunc.c
else
extfunc.so: extfunc.c
	@echo ERROR : cannot compile extfunc - either because gcc not present or gcc command needs to be rewritten for your OS. Check the Makefile.
endif

../../Packages/LanHEP/lhep:
	$(MAKE) -C ../../Packages/LanHEP


clean:
	rm *.so
