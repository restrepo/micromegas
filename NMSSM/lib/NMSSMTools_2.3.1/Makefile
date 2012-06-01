.PHONY: clean all

ifneq ($(MAKECMDGOALS),clean)
include  ../../../CalcHEP_src/FlagsForMake
endif

all : ../nmhdecay  ../nmspec

../nmhdecay:	nmhdecay.F  sources/nmhdecay.a
	$(FC) $(FFLAGS) -o $@ $^ 

../nmspec: nmspec.F sources/nmhdecay.a
	$(FC) $(FFLAGS) -o $@ $^ 

sources/nmhdecay.a: sources/*.f
	$(MAKE) -C sources

clean:
	$(MAKE) -C sources clean
	rm -f  ../nmhdecay ../nmspec 
