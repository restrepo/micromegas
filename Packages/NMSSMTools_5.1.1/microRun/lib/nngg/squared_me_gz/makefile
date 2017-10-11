LIB = ../nngg.a

OBJS = \
  gzabbr1_s.o \
  gzabbr1_angle.o \
  gzabbr1_hel1.o \
  gzabbr1_hel2.o \
  gzvert_neu5.o \
  gzvert_neu5_quu5.o \
  gzvert_neu5_qud5.o \
  gzvert_lpd5_neu5.o \
  gzvert_cha5_neu5.o \
  gzvert.o \
  gzvert_sqv5.o \
  gzvert_squ5.o \
  gzvert_sqe5.o \
  gzvert_sqd5.o \
  gzvert_sle5.o \
  gzvert_sld5.o \
  gzvert_quu5.o \
  gzvert_qud5.o \
  gzvert_lpd5.o \
  gzvert_cha5_cha6.o \
  gzvert_his5.o \
  gzvert_his5_sqv5.o \
  gzvert_his5_squ5.o \
  gzvert_his5_sqe5.o \
  gzvert_his5_sqd5.o \
  gzvert_his5_sle5.o \
  gzvert_his5_sld5.o \
  gzvert_his5_quu5.o \
  gzvert_his5_qud5.o \
  gzvert_his5_lpd5.o \
  gzvert_cha5_cha6_his5.o \
  gzvert_hia5_quu5.o \
  gzvert_hia5_qud5.o \
  gzvert_hia5_lpd5.o \
  gzvert_cha5_cha6_hia5.o \
  gzbox_quu51.o \
  gzbox_quu52.o \
  gzbox_qud51.o \
  gzbox_qud52.o \
  gzbox_lpd51.o \
  gzbox_lpd52.o \
  gzbox_cha51.o \
  gzbox_cha52.o \
  gzbox_cha5_cha61.o \
  gzbox_cha5_cha62.o \
  gzbox_cha5_cha63.o \
  gzbox_cha5_cha64.o \
  gzSquaredME.o

$(LIB): $(LIB)($(OBJS))

$(LIB)($(OBJS)): vars.h renconst.h model.h util.h looptools.h

clean distclean:
	$(RM) $(OBJS) $(LIB)

