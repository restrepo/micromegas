.PHONY: clean

ifneq ($(MAKECMDGOALS),clean)
include ../../../CalcHEP_src/FlagsForMake
endif

cLib = $(CALCHEP)/lib

LIB = ../sources/sources.a 


#ifneq ($(LHAPDFPATH),)
#  LIB += -L$(LHAPDFPATH)  -lLHAPDF $(cLib)/dummy.a
#  DLSET = export LD_RUN_PATH=$(LHAPDFPATH);
#else 
#  LIB += $(cLib)/dummy.a 
#  DLSET= 
#endif   



EXEC = nmhdecay nmspec nmgmsb nmhdecay_grid nmhdecay_rand nmhdecay_mcmc nmhdecayCPV nmhdecayCPV_grid nmhdecayCPV_rand nmhdecayCPV_mcmc nmspec_grid nmspec_rand nmspec_mcmc nmgmsb_grid nmgmsb_rand nmgmsb_mcmc

all: $(EXEC)

%.o: %.f
	$(FC) $(FFLAGS) -o $@ -c $<

nmhdecay:	nmhdecay.o  $(LIB)
	$(FC) $(FFLAGS) -o $@ $^ $(lDL) $(LX11) -lpthread
	
nmhdecay_grid:	nmhdecay_grid.o  $(LIB)
	$(FC) $(FFLAGS) -o $@ $^ $(lDL) $(LX11) -lpthread

nmhdecay_rand:	nmhdecay_rand.o $(LIB)
	$(FC) $(FFLAGS) -o $@ $^ $(lDL) $(LX11) -lpthread

nmhdecay_mcmc:	nmhdecay_mcmc.o $(LIB)
	$(FC) $(FFLAGS) -o $@ $^ $(lDL) $(LX11) -lpthread

nmhdecayCPV:	nmhdecayCPV.o  $(LIB)
	$(FC) $(FFLAGS) -o $@ $^ $(lDL) $(LX11) -lpthread
	
nmhdecayCPV_grid:	nmhdecayCPV_grid.o  $(LIB)
	$(FC) $(FFLAGS) -o $@ $^ $(lDL) $(LX11) -lpthread
	
nmhdecayCPV_rand:	nmhdecayCPV_rand.o  $(LIB)
	$(FC) $(FFLAGS) -o $@ $^ $(lDL) $(LX11) -lpthread
	
nmhdecayCPV_mcmc:	nmhdecayCPV_mcmc.o  $(LIB)
	$(FC) $(FFLAGS) -o $@ $^ $(lDL) $(LX11) -lpthread
	
nmspec: 	nmspec.o $(LIB)
	$(FC) $(FFLAGS) -o $@ $^ $(lDL) $(LX11) -lpthread

nmspec_grid: 	nmspec_grid.o $(LIB)
	$(FC) $(FFLAGS) -o $@ $^ $(lDL) $(LX11) -lpthread

nmspec_rand: 	nmspec_rand.o $(LIB)
	$(FC) $(FFLAGS) -o $@ $^ $(lDL) $(LX11) -lpthread

nmspec_mcmc: 	nmspec_mcmc.o $(LIB)
	$(FC) $(FFLAGS) -o $@ $^ $(lDL) $(LX11) -lpthread

nmgmsb: 	nmgmsb.o $(LIB)
	$(FC) $(FFLAGS) -o $@ $^ $(lDL) $(LX11) -lpthread

nmgmsb_grid: 	nmgmsb_grid.o $(LIB)
	$(FC) $(FFLAGS) -o $@ $^ $(lDL) $(LX11) -lpthread

nmgmsb_rand: 	nmgmsb_rand.o $(LIB)
	$(FC) $(FFLAGS) -o $@ $^ $(lDL) $(LX11) -lpthread

nmgmsb_mcmc: 	nmgmsb_mcmc.o $(LIB)
	$(FC) $(FFLAGS) -o $@ $^ $(lDL) $(LX11) -lpthread

clean:
	rm -f *.o $(EXEC)
