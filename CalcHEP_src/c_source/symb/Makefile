
include ../../FlagsForMake

Idir= chep_crt service2  polynom getmem  ..

override CFLAGS += -Iinclude $(patsubst %,-I../%/include ,$(Idir)) $(SNUM)

special = s_calchep.c V_and_P.c

VPATH = out diagram
OBJ = $(notdir $(subst .c,.o,$(filter-out $(special),$(wildcard *.c diagram/*.c out/*.c))))

bin = ../../bin
lib = ../../lib

all: $(bin)/s_calchep $(bin)/make_VandP 

$(bin)/s_calchep: s_calchep.c $(lib)/symb.a $(lib)/serv.a $(lib)/dynamic_me.a
	$(CC) $(CFLAGS) -o $@ $^  $(LX11) $(lDL) -lm

$(lib)/symb.a:$(lib)/symb.a($(OBJ))
	$(RANLIB) $@
	@rm -f $(OBJ)

$(bin)/make_VandP: make_VandP.c $(lib)/symb.a $(lib)/serv.a
	$(CC) $(CFLAGS) -o $@ $^  $(LX11) -lm
