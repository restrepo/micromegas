LIB = ../nngg.a

OBJS = \
  abbr1_s.o \
  abbr1_angle.o \
  abbr1_hel1.o \
  abbr1_hel2.o \
  vert.o \
  vert_sqv5.o \
  vert_squ5.o \
  vert_sqe5.o \
  vert_sqd5.o \
  vert_sle5.o \
  vert_sld5.o \
  vert_quu5.o \
  vert_qud5.o \
  vert_lpd5.o \
  vert_cha5.o \
  vert_his5.o \
  vert_his5_sqv5.o \
  vert_his5_squ5.o \
  vert_his5_sqe5.o \
  vert_his5_sqd5.o \
  vert_his5_sle5.o \
  vert_his5_sld5.o \
  vert_his5_quu5.o \
  vert_his5_qud5.o \
  vert_his5_lpd5.o \
  vert_cha5_his5.o \
  vert_hia5_quu5.o \
  vert_hia5_qud5.o \
  vert_hia5_lpd5.o \
  vert_cha5_hia5.o \
  box_quu5.o \
  box_qud5.o \
  box_lpd5.o \
  box_cha51.o \
  box_cha52.o \
  SquaredME.o

$(LIB): $(LIB)($(OBJS))

$(LIB)($(OBJS)): vars.h renconst.h model.h util.h looptools.h

clean distclean:
	$(RM) $(OBJS) $(LIB)

