include ../../config.mk

CXXFLAGS+= -I. -I.. -I../../include
LDLIBS=-lboost_filesystem -lboost_iostreams -lboost_system

ifdef BOOST_INCLUDE
    CXXFLAGS+=$(BOOST_INCLUDE)
endif

ifdef BOOST_LIB
    LDFLAGS+=$(BOOST_LIB)
endif

PROGS=elastic-cluster-seq \
      elastic-finalize-seq

all: $(PROGS)

install: $(PROGS)
	mkdir -p ../../bin/
	cp $(PROGS) ../../bin/

clean:
	rm -rf $(PROGS) *.o *.a

.cpp:
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $< -o $@ $(LDLIBS)
