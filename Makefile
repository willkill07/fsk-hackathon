CXX := pgc++
CXXFLAGS := -Minfo=acc -O3 -acc -ta=nvidia,cc50,7.5

SRC := fsk.cpp

.PHONY: all clean

all : fsk

CXX := icpc
CXXFLAGS := -xHOST -O3 -openmp -restrict

clean :
	-rm -rf $(wildcard *.o)
