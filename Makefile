CXX :=
CXXFLAGS :=
LANGTYPE ?= openacc
TARGET ?= gpu
LOADLIBES := -lrt
HOST := $(shell hostname)

ifeq ($(LANG_TYPE),openacc)
# OpenACC Setup
  COMPILER := pgi
  CXXFLAGS += -acc
  SRC := fsk.cpp
  ifeq ($(TARGET),cpu)
    CXXFLAGS += -ta=multicore
  else
    ifeq ($(HOST),liz)
      CXXFLAGS += -ta=nvidia,cc50,7.5
    else
      CXXFLAGS += -ta=nvidia,cc35,7.0
    endif
  endif
else
  SRC := fsk_cpu.cpp
  CXXFLAGS += -std=c++11
endif

ifeq ($(COMPILER),pgi)
# PGI Options
  CXX := pgc++
  CXXFLAGS += -O3 -Minfo=acc,opt,mp
else
# Other Compilers
  ifneq ($(HOST),liz)
    CXX := g++
    CXXFLAGS += -march=native -O3 -fopenmp
  else
    CXX := icpc
    CXXFLAGS += -xHOST -O3 -openmp
  endif
endif

BIN := $(SRC:.cpp=)

.PHONY: all clean cleanall default

default : $(BIN)

all :
	make COMPILER=$(COMPILER) LANG_TYPE=openmp TARGET=cpu
	make COMPILER=pgi LANG_TYPE=openacc TARGET=gpu && mv fsk fsk_acc_gpu
	make COMPILER=pgi LANG_TYPE=openacc TARGET=cpu && mv fsk fsk_acc_cpu

clean :
	-rm -rf fsk_cpu fsk_acc_gpu fsk_acc_cpu
