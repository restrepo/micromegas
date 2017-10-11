
all:HiggsSignals/HiggsSignals HiggsBounds/HiggsBounds

HiggsBounds/HiggsBounds: HiggsBounds/configure
	cd HiggsBounds; ./configure; make

HiggsSignals/HiggsSignals: HiggsBounds/HiggsBounds
	cd HiggsSignals; ./configure; make

HiggsBounds/configure: 
ifneq (, $(shell command -v curl 2> /dev/null))
	curl https://lapth.cnrs.fr/micromegas/downloadarea/packages/hBandS.tgz  -O
else 
ifneq (, $(shell command -v wget 2> /dev/null))
	wget https://lapth.cnrs.fr/micromegas/downloadarea/packages/hBandS.tgz
else 
	$(error "Neither wget nor curl are available, please install wget or curl or change hBandS.make accordingly.")
endif
endif
	tar -xvzf hBandS.tgz 
	rm hBandS.tgz  