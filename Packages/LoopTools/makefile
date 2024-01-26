# --- variables defined by configure ---

SRC = ./src
QUADSUFFIX = 
PREFIX = x86_64-Linux
LIBDIRSUFFIX = 64

EXE = 
DEF = -D
NOUNDERSCORE = 0
QUAD = 0
QUADSIZE = 16

FC = gfortran
FFLAGS = -O3 -ffixed-line-length-none -fno-range-check \
  $(DEF)QUAD=$(QUAD) $(DEF)QUADSIZE=$(QUADSIZE) \
  $(DEF)U77EXT=0

CDEFS = -DQUAD=$(QUAD) -DQUADSIZE=$(QUADSIZE) \
  -DNOUNDERSCORE=$(NOUNDERSCORE) \
  -DBIGENDIAN=0

CC = gcc
CFLAGS = -O3 -g -fomit-frame-pointer -ffast-math -fPIC -Wall -m64

CXX = g++
CXXFLAGS = -O3 -g -fomit-frame-pointer -ffast-math -fPIC -Wall -m64

ML = 1
MCC = mcc
MCDEFS = -D__int64=long\ long\ int
MCFLAGS = -n
MCLIBS = -lpthread -lrt

LDFLAGS = -L/usr/libexec/gcc/x86_64-redhat-linux/8/liblto_plugin.so -L/usr/lib/gcc/x86_64-redhat-linux/8 -L/usr/lib/gcc/x86_64-redhat-linux/8/../../../../lib64 -L/lib/../lib64 -L/usr/lib/../lib64 -L/usr/lib/gcc/x86_64-redhat-linux/8/../../.. -lgfortran -lm -lgcc_s -lgcc -lquadmath -lm -lgcc_s -lgcc -lgcc_s -lgcc

PARALLEL = -j 4

# --- end defs by configure ---

BLD = build$(QUADSUFFIX)

LIBDIR = $(PREFIX)/lib$(LIBDIRSUFFIX)
INCLUDEDIR = $(PREFIX)/include
BINDIR = $(PREFIX)/bin

LIB = libooptools$(QUADSUFFIX).a
FE = lt$(QUADSUFFIX)$(EXE)
MFE = LoopTools$(QUADSUFFIX)$(EXE)
FCC = fcc$(QUADSUFFIX)
FXX = f++$(QUADSUFFIX)
INCLUDE = $(BLD)/looptools.h $(BLD)/clooptools.h

ARGS = $(PARALLEL) \
  LIB="$(LIB)" \
  FE="$(FE)" \
  MFE="$(MFE)" \
  FCC="$(FCC)" \
  FXX="$(FXX)" \
  EXE="$(EXE)" \
  DEF="$(DEF)" \
  NOUNDERSCORE="$(NOUNDERSCORE)" \
  XFC="$(FC) $(FFLAGS) -I." \
  F90="$(F90)" \
  CDEFS="$(CDEFS)" \
  MCDEFS="$(MCDEFS)" \
  CC="$(CC)" \
  CFLAGS="$(CFLAGS)" \
  CXX="$(CXX)" \
  CXXFLAGS="$(CXXFLAGS)" \
  ML="$(ML)" \
  MCC="$(MCC)" \
  MCFLAGS="$(MCFLAGS)" \
  MCLIBS="$(MCLIBS)" \
  AR="$(AR)" \
  RANLIB="$(RANLIB)" \
  NM="$(NM)" \
  DLLTOOL="$(DLLTOOL)" \
  LDFLAGS="$(LDFLAGS)" \
  LIBPATH="$(LIBPATH)"


default all lib frontend mma: force
	cd $(BLD) && $(MAKE) $(ARGS) $@

install: lib frontend
	-mkdir $(PREFIX)
	-mkdir $(LIBDIR) $(BINDIR) $(INCLUDEDIR)
	cp -p $(BLD)/$(LIB) $(LIBDIR)
	cp -p $(INCLUDE) $(INCLUDEDIR)
	strip $(BLD)/$(FE)
	cp -p $(BLD)/$(FCC) $(BLD)/$(FE) $(BINDIR)
	rm -f $(BINDIR)/$(FXX)
	ln -s $(FCC) $(BINDIR)/$(FXX)
	test ! -f $(BLD)/$(MFE) || { \
	  test -z "$(EXE)" || ldd $(BLD)/$(MFE) | awk '/\/usr\/bin\/cyg/ { system("cp -p " $$3 " $(BINDIR)/"); }' ; \
	  strip $(BLD)/$(MFE) ; \
	  cp -p $(BLD)/$(MFE) $(BINDIR); \
	}

force: $(BLD)/timestamp

$(BLD)/timestamp:
	-mkdir $(BLD)
	find "`cd $(SRC) ; pwd`" -follow -exec ln -sf {} $(BLD) \;
	touch $(BLD)/timestamp

clean:
	rm -fr $(BLD)

