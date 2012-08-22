include ../../config.mk

CXXFLAGS=-O3 -I. -I.. -I../../include
LDFLAGS=-L../
LDLIBS=-liomanip $(MPEFLAGS)

ifdef BOOST_INCLUDE
    CXXFLAGS+=$(BOOST_INCLUDE)
endif

PROGS=elastic-sketch-mpi

all: $(PROGS)

install: $(PROGS)
	mkdir -p ../../bin/
	cp $(PROGS) ../../bin/

clean:
	rm -rf $(PROGS)

.cpp:
	$(MPICXX) $(CXXFLAGS) $(LDFLAGS) $< -o $@ $(LDLIBS)
