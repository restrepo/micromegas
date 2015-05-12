.PHONY: clean

ifneq ($(MAKECMDGOALS),clean)
include ../../../CalcHEP_src/FlagsForMake
endif


LIB =  ../sources/nmhdecay.a 

ALIB =  ../sources/nmhdecay.a $  

EXEC = nmhdecay nmspec 

all: $(EXEC)

%.o: %.f
	$(FC) $(FFLAGS) -o $@ -c $<

nmhdecay:	nmhdecay.o  $(LIB)
	$(FC) $(FFLAGS) -o $@ $^ $(lDL) $(LX11)
	
nmspec: 	nmspec.o $(LIB)
	$(FC) $(FFLAGS) -o $@ $^ $(lDL) $(LX11)

clean:
	rm -f *.o $(EXEC)
