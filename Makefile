all:	basic omp mpi

install:	omp-install mpi-install

basic:
	$(MAKE) -C src/

omp:
	$(MAKE) -C src/omp

mpi: basic
	$(MAKE) -C src/mpi

omp-install: omp
	$(MAKE) -C src/omp install

mpi-doc: mpi-install
	mkdir -p doc/
	bin/elastic-sketch-mpi > doc/elastic-sketch-mpi.txt

mpi-install: mpi
	$(MAKE) -C src/mpi install


clean:
	$(MAKE) -C src/ clean
	$(MAKE) -C src/omp clean
	$(MAKE) -C src/mpi clean

distclean: clean
	rm -rf bin/
