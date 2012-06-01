.PHONY: all clean

# read compiler flags
ifneq ($(MAKECMDGOALS),clean)
include ../../CalcHEP_src/FlagsForMake
endif

# path to your project files.


# all: mLib.$(SO) aLib.a 

# 

all: aLib.a

aLib.a:aLib.a(rdVarlHiggs.o)
	$(RANLIB) aLib.a
	rm -f $(rdVarlHiggs.o)

clean: 
	rm -f *.o *.a *.so *.dll  *.exe so_located
