include ../../config.mk

CXXFLAGS+= -I. -I.. -I../../include
LDLIBS=-lboost_filesystem -lboost_iostreams -lboost_system

ifdef BOOST_INCLUDE
    CXXFLAGS+=$(BOOST_INCLUDE)
endif

ifdef BOOST_LIB
    LDFLAGS+=$(BOOST_LIB)
endif

PROGS=elastic-prepare

all: $(PROGS)
	mkdir -p ../../bin/
	cp $(PROGS) ../../bin/

install: $(PROGS)
	mkdir -p $(PREFIX)/bin
	cp $(PROGS) $(PREFIX)/bin

clean:
	rm -rf $(PROGS) *.o *.a

.cpp:
	$(CXX) $(OMPFLAGS) $(CXXFLAGS) $(LDFLAGS) $< -o $@ $(LDLIBS)
