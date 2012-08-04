include ../config.mk

CXXFLAGS=-O3

libs: iomanip.o
	$(AR) rcs libiomanip.a iomanip.o

clean:
	rm -rf *.o *.a

.cpp.o:
	$(CXX) -c $(CXXFLAGS) $< -o $@
