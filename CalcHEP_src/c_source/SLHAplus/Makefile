.PHONY: all

include ../../FlagsForMake

all: $(CALCHEP)/lib/libSLHAplus.a


OBJ =  $(notdir $(subst .c,.o,$(wildcard *.c))) faux.o
override CFLAGS += -I./include

$(CALCHEP)/lib/libSLHAplus.a:$(CALCHEP)/lib/libSLHAplus.a($(OBJ))
	$(RANLIB) $(CALCHEP)/lib/libSLHAplus.a
	rm -f *.o
 
