# --- variables defined by configure ---

ifneq ($(MAKECMDGOALS),clean)
include ../../CalcHEP_src/FlagsForMake
endif 

FFLAGS +=   -ffixed-line-length-none

.PHONY: force squared_me_gg squared_me_gz  renconst util

export VPATH := $(CURDIR):$(CURDIR)/renconst

INCLUDE := $(patsubst %,-I%,$(subst :, ,$(VPATH)))

FFLAGS += $(INCLUDE)

ifdef DEBUG
FFLAGS += -g -DDEBUG
endif

export FC
export FFLAGS


lGamma.exe : main.f nngg.a LoopTools-2.1/libooptools.a 
	$(FC) $(FFLAGS) -o lGamma.exe  main.f  nngg.a LoopTools-2.1/libooptools.a ../../CalcHEP_src/lib/libSLHAplus.a 

LoopTools-2.1/libooptools.a:
	$(MAKE) -C LoopTools-2.1
#subdir = squared_me_gg squared_me_gz renconst  util 
subdir = util
#OBJ= $(notdir $(subst .c,.o, $(filter-out $(special),$(wildcard *.c diagram/*.c out/*.c))))

nnggOBJ = $(subst .F,.o, $(wildcard *.F squared_me_gg/*.F  squared_me_gz/*.F  renconst/*.F))
#nngg.a : $(subdir) nngg.a(run_gg.o run_gz.o mdl_ini.o interface.o interfacegz.o)

nngg.a : $(subdir) nngg.a($(nnggOBJ))
nngg.a(run.o): 2to2.h renconst.h looptools.h model.h

nngg.a(mdl_ini.o):model.h

$(subdir) : 
	$(MAKE) -C  $@

clean:	
	$(MAKE) -C LoopTools-2.1 clean
	rm -f nngg.a lGamma.exe *.o */*.o
