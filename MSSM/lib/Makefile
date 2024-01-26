
SUSPECT= $(CURDIR)/../../Packages/SuSpect_2.41

#read compiler options ..

ifneq ($(MAKECMDGOALS),clean)
include ../../CalcHEP_src/FlagsForMake
endif 

.PHONY: all clean

VPATH = SuSpect:../../sources

SuSpectSRC = $(wildcard SuSpect/*.f)


aLibOBJ = $(subst .c,.o,$(filter-out $(special),$(wildcard *.c))) 

ifneq ($(FC),)
   all:  suspect.exe  aLib.a 
   special = isajetLesH.c
else
   all: aLib.a
   special = isajetLesH.c suspect_call.c
endif

suspect.exe: $(SuSpectSRC) 
	$(MAKE) -C $(SUSPECT)
	echo \#define SUSPECT \"$(SUSPECT)\"  > suspect_path.h
aLib.a:aLib.a($(aLibOBJ))
	$(RANLIB) aLib.a
	rm -f $(aLibOBJ)

# for Isajet before 7.79 add  flag  -Dlt779 for FC
isajet_slha: isasugra_slha.F isajetLesH.c
	$(CC) $(CFLAGS)   -c  isajetLesH.c
	$(FC) $(FFLAGS) -o isajet_slha  isasugra_slha.F isajetLesH.o \
   -L$(ISAJET) -L$(CERNLIB) -lisajet -lmathlib 

clean: 
	rm -f *.o *.a  *.so *.dll *.exe *_path.h isajet_slha so_locations 
	$(MAKE) -C nngg12     clean