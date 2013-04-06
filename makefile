# UjoImro, 2013

# ## absolute directory for the subprojects
export BASEDIR=$(shell pwd)

export CXXFLAGS=-std=c++0x -I$(BASEDIR)/opencl
export LDFLAGS=-lstdc++ -lOpenCL -lm

# optimization flags
CXXFLAGS+=-O0 -g 

# includes 
CXXFLAGS+=-I/home/ujoimro/Inst/opencv/build/OpenCV-2.4.2/optimized/install/include

CFLAGS+=-Wimplicit -Werror=implicit-function-declaration -Werror=int-to-pointer-cast

# c++0x and c++ specific flags
CXXFLAGS+=-std=c++0x

## Memory Debug Mode
#LDFLAGS=-lboost_serialization -L/home/ujoimro/Inst/opencv/build/OpenCV-2.4.2/optimized/install/lib -lopencv_core -lstdc++ -lm -lefenc
LDFLAGS+=-lboost_serialization -L/home/ujoimro/Inst/opencv/build/OpenCV-2.4.2/optimized/install/lib -lopencv_core 


all:
	+make -C opencl
	+make -C base
	+make -C mlp

clean:
	make -C opencl clean
	make -C base clean
	make -C mlp clean

# LuM end of file
