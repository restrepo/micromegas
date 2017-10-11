
.PHONY: all

all: *.lhp ../../Packages/LanHEP/lhep
	../../Packages/LanHEP/lhep -ca -evl 2 idublet.lhp

../../Packages/LanHEP/lhep:
	$(MAKE) -C ../../Packages/LanHEP
	