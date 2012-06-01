.PHONY: all clean

ifneq ($(MAKECMDGOALS),clean)
include ../../CalcHEP_src/FlagsForMake 
endif



all:  aLib.a

mLib.$(SO):lzp.c bessj0.c bessj1.c bessy0.c bessy1.c bessjy.c chebev.c nrutil.c simp.c couplings.c dd.c
	$(CC) $(CFLAGS)  $(SHARED) -o $(target) $^ -lm

aLib.a:aLib.a(rdVarRHNM.o lzp.o bessj0.o bessj1.o bessy0.o bessy1.o bessjy.o chebev.o nrutil.o simp.o couplings.o dd.o)
	$(RANLIB) aLib.a
	rm -f rdVarRHNM.o

clean: 
	rm -f *.o *.a 
