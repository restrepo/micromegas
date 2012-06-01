
.PHONY: all libs clean

# read compiler flags
ifneq ($(MAKECMDGOALS),clean)
AllFlags = ../CalcHEP_src/FlagsForMake
ifeq (,$(wildcard $(AllFlags) )) 
$(error File $(AllFlags) is absent. Presumably you forgot to compile main code)
endif 
include ../CalcHEP_src/FlagsForMake 
endif 

cLib = $(CALCHEP)/lib

# files to compile
SSS =  ../sources/micromegas.a  $(cLib)/dynamic_me.a ../sources/micromegas.a   work/work_aux.a \
  $(wildcard lib/aLib.a)  $(cLib)/sqme_aux.$(SO) $(cLib)/libSLHAplus.a   \
 ../CalcHEP_src/lib/num_c.a  ../CalcHEP_src/lib/serv.a $(LX11)  

all: $(main) libs
ifeq (,$(main)) 
	@echo Main program is not specified. Use gmake main='<code of main program>'
else  
  ifeq (.cpp,$(suffix $(main)))
	$(CXX) $(CXXFLAGS)  -o $(main:.cpp=) $(main) $(SSS)  $(lDL)  -lm
  else  
    ifeq (.c,$(suffix $(main)))
	$(CC)  $(CFLAGS) -o $(main:.c=) $(main)  $(SSS)  $(lDL)  -lm 
    else  
      ifeq (.F,$(suffix $(main)))
        ifeq ($(FC), )
	$(error  Fortran compiler was not detected.)
        else 
	$(FC) $(FFLAGS) -o $(main:.F=) $(main)  $(SSS)  $(lDL)  -lm
        endif
      else
	$(error Only .c, .cpp, and .F  main programs  are supported).
      endif
    endif
  endif  
endif



libs:
	$(MAKE) -C work
	$(MAKE) -C lib

clean: 
	$(MAKE) -C lib  clean
	$(MAKE) -C work clean 
	rm -f calchep/tmp/* calchep/results/*
	../sources/cleanexe
	rm -f LesHin LesHout suspect2_lha.in suspect2_lha.out suspect2.out tab
	rm -f Control.in LesHouches.in Messages.out SPheno.out slha.in slha.out
