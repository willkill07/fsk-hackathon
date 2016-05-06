CXX :=
CXXFLAGS :=
LANGTYPE ?= openacc
TARGET ?= gpu
LOADLIBES := -lrt
HOST := $(shell hostname)

#TAU := tau

ifeq ($(LANG_TYPE),openacc)
# OpenACC Setup
  COMPILER := pgi
  CXXFLAGS += -acc
  SRC := fsk.cpp
  ifeq ($(TARGET),cpu)
    CXXFLAGS += -ta=multicore
  else
    ifeq ($(HOST),liz)
      CXXFLAGS += -ta=nvidia:managed,maxwell,cuda7.5,keepptx,ptxinfo,maxregcount:32
    else
      CXXFLAGS += -ta=nvidia,kepler,cuda7.0,keepptx,ptxinfo,maxregcount:32
    endif
  endif
else
  SRC := fsk_cpu.cpp
  CXXFLAGS += -std=c++11
endif

ifeq ($(COMPILER),pgi)
# PGI Options
  CXX := $(TAU) pgc++
  CXXFLAGS += -O3 -Minfo=acc,opt,mp
else
# Other Compilers
  ifneq ($(HOST),liz)
    CXX := g++
    CXXFLAGS += -march=native -O3 -fopenmp
  else
    CXX := $(TAU) icpc
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
