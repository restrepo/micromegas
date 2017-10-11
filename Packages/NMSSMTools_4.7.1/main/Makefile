.PHONY: clean

ifneq ($(MAKECMDGOALS),clean)
include ../../../CalcHEP_src/FlagsForMake
endif

LIB = ../sources/nmhdecay.a 

ALIB = ../sources/nmhdecay.a 


EXEC = nmhdecay nmspec

all: $(EXEC)

%.o: %.f
	$(FC) $(FFLAGS) -o $@ -c $<

nmhdecay:	nmhdecay.o  $(LIB)
	$(FC) $(FFLAGS) -o $@ $^ $(lDL) $(LX11)
	
nmhdecay_grid:	nmhdecay_grid.o  $(LIB)
	$(FC) $(FFLAGS) -o $@ $^ $(lDL) $(LX11)

nmhdecay_rand:	nmhdecay_rand.o $(LIB)
	$(FC) $(FFLAGS) -o $@ $^ $(lDL) $(LX11)

nmhdecay_mcmc:	nmhdecay_mcmc.o $(LIB)
	$(FC) $(FFLAGS) -o $@ $^ $(lDL) $(LX11)

nmspec: 	nmspec.o $(LIB)
	$(FC) $(FFLAGS) -o $@ $^ $(lDL) $(LX11)

nmspec_grid: 	nmspec_grid.o $(LIB)
	$(FC) $(FFLAGS) -o $@ $^ $(lDL) $(LX11)

nmspec_rand: 	nmspec_rand.o $(LIB)
	$(FC) $(FFLAGS) -o $@ $^ $(lDL) $(LX11)

nmspec_mcmc: 	nmspec_mcmc.o $(LIB)
	$(FC) $(FFLAGS) -o $@ $^ $(lDL) $(LX11)

nmgmsb: 	nmgmsb.o $(LIB)
	$(FC) $(FFLAGS) -o $@ $^ $(lDL) $(LX11)

nmgmsb_grid: 	nmgmsb_grid.o $(LIB)
	$(FC) $(FFLAGS) -o $@ $^ $(lDL) $(LX11)

nmgmsb_rand: 	nmgmsb_rand.o $(LIB)
	$(FC) $(FFLAGS) -o $@ $^ $(lDL) $(LX11)

nmgmsb_mcmc: 	nmgmsb_mcmc.o $(LIB)
	$(FC) $(FFLAGS) -o $@ $^ $(lDL) $(LX11)

clean:
	rm -f *.o $(EXEC)
