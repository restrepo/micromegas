LIB = ../nngg.a

OBJS = \
  abbr1s.o \
  abbr1angle.o \
  abbr1hel.o \
  vert.o \
  box.o \
  SquaredME.o

$(LIB): $(LIB)($(OBJS))

$(LIB)($(OBJS)): vars.h 

#renconst.h model.h util.h looptools.h

clean distclean:
	$(RM) $(OBJS) $(LIB)

