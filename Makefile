CXX := pgc++
CXXFLAGS := -Minfo=acc -O3 -acc -ta=nvidia,cc35,7.0

SRC := fsk.cpp

.PHONY: all clean

all : fsk

clean :
	-rm -rf $(wildcard *.o)
