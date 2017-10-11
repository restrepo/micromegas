
.PHONY:  clean all

ifneq ($(MAKECMDGOALS),clean)
include ../../CalcHEP_src/FlagsForMake
endif

all: extfunc work_aux.a 

platform=$(shell uname)
workOBJ =  VandP.o  path.o  extfunc.o
override CFLAGS += -I../../CalcHEP_src/include

ifeq ($(platform),Darwin)
extfunc: extfunc.c
	gcc -c extfunc.c
	ar rcs extfunc.a extfunc.o
	sed -e "s/extfunc.*/extfunc.a/g" models/extlib1.mdl > models/extlib1.mdl_new && mv models/extlib1.mdl_new models/extlib1.mdl
	cp models/extlib1.mdl ../calchep/models/extlib1.mdl
else ifeq ($(platform),Linux)
extfunc: extfunc.c
	gcc -fPIC -shared -o extfunc.so extfunc.c
	sed -i -e "s/extfunc.*/extfunc.so/g" models/extlib1.mdl ../calchep/models/extlib1.mdl
else
extfunc: extfunc.c
	gcc -fPIC -shared -o extfunc.so extfunc.c
	sed -i -e "s/extfunc.*/extfunc.so/g" models/extlib1.mdl ../calchep/models/extlib1.mdl
	@echo If gcc/sed commands are not working : ERROR, cannot compile extfunc - either because gcc/sed not present or gcc/sed commands need to be rewritten for your OS. Check the Makefile.
endif

work_aux.a: work_aux.a($(workOBJ))
	$(RANLIB) work_aux.a 
	rm -f $(workOBJ)

VandP.c: models/vars1.mdl models/func1.mdl models/prtcls1.mdl
	../../CalcHEP_src/bin/make_VandP models 1 
	@echo
	@echo - - - - - Please wait, the following step may take a few seconds - - - - - 

path.c:
	echo char \* WORK=\"`pwd`\"\; >  path.c


clean:
	rm -f *.a *.o *.so lock_ LOCK  tmp/*  so_generated/* lock_
	rm -f  path.c VandP.c autoprot.h
	rm -rf _*_*_ results/*
	sed -e "s/extfunc.*/extfunc/g" models/extlib1.mdl > models/extlib1.mdl_new && mv models/extlib1.mdl_new models/extlib1.mdl
	cp models/extlib1.mdl ../calchep/models/extlib1.mdl