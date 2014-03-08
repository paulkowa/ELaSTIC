include config.mk

all:	basic seq omp mpi

install:	seq-install omp-install mpi-install
	mkdir -p $(PREFIX)/share/doc/elastic-`cat VERSION`
	cp doc/* $(PREFIX)/share/doc/elastic-`cat VERSION`

basic:
	$(MAKE) -C src/

seq: basic
	$(MAKE) -C src/seq

omp: basic
	$(MAKE) -C src/omp

mpi: basic
	$(MAKE) -C src/mpi

seq-install: seq
	$(MAKE) -C src/seq install

omp-install: omp
	$(MAKE) -C src/omp install

mpi-doc: mpi-install
	mkdir -p doc/
	bin/elastic-sketch > doc/elastic-sketch.txt

mpi-install: mpi
	$(MAKE) -C src/mpi install

others-install: seq-install omp-install

sketch: mpi
sketch-doc: mpi-doc
sketch-install: mpi-install

clean:
	$(MAKE) -C src/ clean
	$(MAKE) -C src/seq clean
	$(MAKE) -C src/omp clean
	$(MAKE) -C src/mpi clean

distclean: clean
	rm -rf bin/

release: all mpi-doc
	$(MAKE) distclean
