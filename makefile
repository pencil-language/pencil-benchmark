# UjoImro, 2013

# ## absolute directory for the subprojects
BASEDIR=$(shell pwd)

CXX=g++
CC=gcc
LD=ld
CXXFLAGS=-std=c++0x -I$(BASEDIR)/opencl -I$(BASEDIR)/core
CFLAGS=-lOpenCL -lm
LDFLAGS=-lstdc++ -lOpenCL -lm -lboost_serialization -L/home/ujoimro/Inst/opencv/build/OpenCV-2.4.2/optimized/install/lib -lopencv_core -lstdc++ -lm

# optimization flags
CXXFLAGS+=-O3 -mtune=native

# includes 
CXXFLAGS+=-I/home/ujoimro/Inst/opencv/build/OpenCV-2.4.2/optimized/install/include

CFLAGS+=-Wimplicit -Werror=implicit-function-declaration -Werror=int-to-pointer-cast -I$(BASEDIR)/opencl -I$(BASEDIR)/core

# c++0x and c++ specific flags
CXXFLAGS+=-std=c++0x

## Memory Debug Mode
#LDFLAGS=-lboost_serialization -L/home/ujoimro/Inst/opencv/build/OpenCV-2.4.2/optimized/install/lib -lopencv_core -lstdc++ -lm -lefenc
# LDFLAGS+=-lboost_serialization -L/opt/local/lib64 -lopencv_core 

export CXX
export CC
export LD
export CFLAGS
export BASEDIR
export LDFLAGS
export CXXFLAGS
export CFLAGS

all:
	+make -C core
	+make -C opencl
	+make -C base
	+make -C mlp

clean:
	make -C core clean
	make -C opencl clean
	make -C base clean
	make -C mlp clean

# LuM end of file
