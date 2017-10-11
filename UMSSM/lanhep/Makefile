#	This Makefile and the UMSSM model files which it compiles are used in LanHEP version 3.2.0
FORTLIB = -lgfortran
SHLIBS = extfunc.so

.PHONY: all clean distclean

all: $(SHLIBS)
	../lhep startup.src -ca -evl 2
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
clean:
	rm *.so
