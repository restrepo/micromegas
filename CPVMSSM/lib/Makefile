.PHONY: all

ifneq ($(MAKECMDGOALS),clean)
include ../../CalcHEP_src/FlagsForMake
endif

all:  CPsuperH2/cpsuperh.exe aLib.a 


OBJ = $(subst .c,.o,$(wildcard *.c)) 



aLib.a:  aLib.a($(OBJ))
	$(RANLIB) aLib.a
	rm -f $(OBJ)

aLib.a(cpsuperh_call.o):localpath.h cpsuperh_call.c

localpath.h:
	echo \#define lPath   \"$(CURDIR)\"  >  localpath.h

	
#edm.o : edm.c
#	 cc -c edm.c
        
CPsuperH2/cpsuperh.exe: CPsuperH2/*.f  HiggsBound/libHB.a
	$(MAKE) -C CPsuperH2

HiggsBound/libHB.a: HiggsBound/*.f HiggsBound/*.F
	$(MAKE) -C HiggsBound

clean: 
	rm -f *.o *.a  localpath.h
	$(MAKE) -C CPsuperH2 clean
	$(MAKE) -C HiggsBound clean
  
