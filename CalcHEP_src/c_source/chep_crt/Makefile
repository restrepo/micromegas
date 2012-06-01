
include ../../FlagsForMake

Idir= getmem

override CFLAGS += $(LOCK) -Iinclude $(patsubst %,-I../%/include ,$(Idir)) $(HX11)

OBJ= $(subst .c,.o,$(wildcard *.c)) X11_crt0.o

ifeq (,$(LX11))
 vpath X11_crt0.c noxwin
else 
 vpath X11_crt0.c xwin
endif

lib = ../../lib
$(lib)/serv.a:$(lib)/serv.a($(OBJ))
	$(RANLIB) $@
	@rm -f $(OBJ)
