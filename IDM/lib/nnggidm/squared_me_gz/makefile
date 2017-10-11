PREFIX = AZ

FFLAGS += -DPREFIX=$(PREFIX)

LIB = ../nngg.a

OBJS = \
  gzabbr1s.o \
  gzabbr1a.o \
  gzabbr1h.o \
  gzvert.o \
  gzbox.o \
  gzSquaredME.o

$(LIB): $(LIB)($(OBJS))

$(LIB)($(OBJS)): vars.h renconst.h model.h util.h looptools.h

clean distclean:
	$(RM) $(OBJS) $(LIB)

