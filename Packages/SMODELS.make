
all:smodels-v1.1.1

smodels-v1.1.1: smodels-v1.1.1.tgz
	tar -xvzf smodels-v1.1.1.tgz
	rm smodels-v1.1.1.tgz
	make -C smodels-v1.1.1

smodels-v1.1.1.tgz:
ifneq (, $(shell command -v curl 2> /dev/null))
#	curl http://theory.sinp.msu.ru/~pukhov/smodels-v1.1.1.tgz  -O
	curl https://lapth.cnrs.fr/micromegas/downloadarea/packages/smodels-v1.1.1.tgz  -O
else 
ifneq (, $(shell command -v wget 2> /dev/null))
#	wget http://theory.sinp.msu.ru/~pukhov/smodels-v1.1.1.tgz
	wget https://lapth.cnrs.fr/micromegas/downloadarea/packages/smodels-v1.1.1.tgz 
else 
	$(error "Neither wget nor curl are available, please install wget or curl or change SMODELS.make accordingly.")
endif
endif
