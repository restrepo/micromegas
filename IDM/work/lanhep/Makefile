
.PHONY: all

all: *.lhp ../../../Packages/LanHEP/lhep
	../../../Packages/LanHEP/lhep -fa -evl 2 idublet.lhp

../../../Packages/LanHEP/lhep:
	$(MAKE) -C ../../../Packages/LanHEP
	