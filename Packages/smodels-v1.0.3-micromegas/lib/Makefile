all: nllfast
clean: clean_nllfast

nllfast: nllfast-1.2 nllfast-2.1
clean_nllfast: clean_nllfast-1.2 clean_nllfast-2.1

nllfast-1.2:
	cd nllfast/nllfast-1.2 && make
clean_nllfast-1.2:
	cd nllfast/nllfast-1.2 && make clean
nllfast-2.1:
	cd nllfast/nllfast-2.1 && make
clean_nllfast-2.1:
	cd nllfast/nllfast-2.1 && make clean

.PHONY:
