
include ../../FlagsForMake

override CFLAGS += -Iinclude -DMAKE=\"$(MAKE)\" -Wmissing-prototypes

special =

#OBJ= $(notdir $(subst .c,.o, $(filter-out $(special),$(wildcard *.c diagram/*.c out/*.c))))

OBJ = dynamic_cs.o kin4.o
lib = ../../lib

all: $(lib)/dynamic_me.a

$(bin)/s_calchep: s_calchep.c $(lib)/symb.a $(lib)/serv.a
	$(CC) $(CFLAGS) -o $@ $^  $(LX11) -lm

$(lib)/dynamic_me.a:$(lib)/dynamic_me.a($(OBJ))
	$(RANLIB) $@
	@rm -f $(OBJ)

