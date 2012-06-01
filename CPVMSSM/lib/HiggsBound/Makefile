ifneq ($(MAKECMDGOALS),clean)
include ../../../CalcHEP_src/FlagsForMake
endif

override FFLAGS += -ffixed-line-length-none -DINSTDIR=\"$(CURDIR)\"

OBJ = utils.o f77subroutines.o cpshInterface.o
libHB.a : libHB.a($(OBJ))
	$(RANLIB) libHB.a
	@rm -f $(OBJ)
	
clean: 
	@rm -f 	$(OBJ) libHB.a

