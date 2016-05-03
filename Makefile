#ifdef CPU
CXX := icpc
CXXFLAGS := -xHOST -O3 -openmp -restrict -std=c++11
SRC := fsk_cpu.cpp
# else
# CXX := pgc++
# CXXFLAGS := -Minfo=acc -O3 -acc -ta=nvidia,cc50,7.5
# SRC := fsk.cpp
# endif

TARGET := $(SRC:.cpp=)

.PHONY: all clean

all : $(TARGET)

clean :
	-rm -rf $(TARGET)
