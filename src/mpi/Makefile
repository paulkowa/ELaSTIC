include ../../config.mk

CXXFLAGS+=-I. -I.. -I../../include
LDFLAGS=-L../
LDLIBS=-liomanip $(MPE_FLAGS)

ifdef BOOST_INCLUDE
    CXXFLAGS+=$(BOOST_INCLUDE)
endif

ifdef MPE_LIB
    LDFLAGS+=$(MPE_LIB)
endif

ifdef WITH_MPE
    CXXFLAGS+=$(MPE_INCLUDE) $(WITH_MPE)
endif

PROGS=elastic-sketch

all: $(PROGS)

install: $(PROGS)
	mkdir -p ../../bin/
	cp $(PROGS) ../../bin/

clean:
	rm -rf $(PROGS) *.o *.a

.cpp:
	$(MPICXX) $(CXXFLAGS) $(LDFLAGS) $< -o $@ $(LDLIBS)
