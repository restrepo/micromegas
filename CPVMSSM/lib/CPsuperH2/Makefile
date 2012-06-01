ifneq ($(MAKECMDGOALS),clean)
include ../../../CalcHEP_src/FlagsForMake
endif

cpsuperh.exe: libcpsuperh.a cpsuperh2.f ../HiggsBound/libHB.a
	$(FC) $(FFLAGS) -o cpsuperh.exe cpsuperh2.f libcpsuperh.a ../HiggsBound/libHB.a

OBJ = $(notdir $(subst .f,.o,$(filter-out cpsuperh2.f,$(wildcard *.f))))
libcpsuperh.a:$(OBJ)

libcpsuperh.a:libcpsuperh.a($(OBJ))
	$(RANLIB) libcpsuperh.a
	rm *.o

clean:
	 rm -f cpsuperh.exe libcpsuperh.a
