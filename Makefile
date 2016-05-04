ifeq ($(LANG_TYPE),openacc)
  CXX := pgc++
  CXXFLAGS := -Minfo=acc -O3 -acc
  SRC := fsk.cpp
  ifeq ($(TARGET),cpu)
    CXXFLAGS += -ta=multicore
  else
    CXXFLAGS += -ta=nvidia,cc35,7.0
  endif
else
  CXX := g++
  CXXFLAGS := -march=native -O3 -fopenmp -std=c++11
  SRC := fsk_cpu.cpp
endif

BIN := $(SRC:.cpp=)

.PHONY: all clean cleanall default

default : $(BIN)

all :
	make LANG_TYPE=openmp TARGET=cpu
	make LANG_TYPE=openacc TARGET=gpu && mv fsk fsk_acc_gpu
	make LANG_TYPE=openacc TARGET=cpu && mv fsk fsk_acc_cpu

clean :
	-rm -rf fsk_cpu fsk_acc_gpu fsk_acc_cpu
