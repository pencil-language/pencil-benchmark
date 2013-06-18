# UjoImro, 2013

# ## absolute directory for the subprojects
BASEDIR=$(shell pwd)

# Use icc by default, but allow the user to it
ifeq ($(USE_GCC), 1)
	CXX=g++
	CC=gcc
	LD=ld
	CXXFLAGS=-std=c++11 -I$(BASEDIR)/opencl -I$(BASEDIR)/core
	LDLIBS+=-lstdc++ -lm
else
	CXX?=icc
	CC?=icc
	LD?=xiar
	CXXFLAGS=-std=c++0x -I$(BASEDIR)/opencl -I$(BASEDIR)/core
	LDLIBS=-lstdc++ -lm
endif


LDFLAGS+=-L/home/ujoimro/Inst/opencv/build/opencv-2.4.5/release/install/lib

LDLIBS+=-lboost_serialization -lOpenCL
LDLIBS+=-lopencv_core -lopencv_ocl -lopencv_imgproc -lopencv_flann -lopencv_highgui -lopencv_features2d -lopencv_objdetect -lopencv_video

ifneq ($(NO_INTEL_LIBS), 1)
	CXXFLAGS+=-I/opt/intel/composer_xe_2013.4.183/ipp/include
	LDLIBS+=-ltbb -lipps -lippi -lmkl_core -lmkl_vml_mc3 -lmkl_vml_avx -lmkl_vml_def -lmkl_vml_mc2 -lmkl_vml_p4n -lmkl_vml_mc -lmkl_vml_cmpt -lmkl_vml_avx2 -lmkl_intel_lp64 -lmkl_sequential
endif

# optimization flags
CXXFLAGS+=-O3

ifeq (CXX, icc)
	CXXFLAGS+=-xHOST -ipo
	CXXFLAGS+=-gcc-name=gcc-4.6
	CXXFLAGS+=-I/opt/intel/composer_xe_2013/ipp/include
endif

ifeq (CXX, g++)
	CXXFLAGS+=-mtune=native
endif

# includes 
CXXFLAGS+=-I/home/ujoimro/Inst/opencv/build/opencv-2.4.5/release/install/include
#CXXFLAGS+=-I/home/ujoimro/Inst/opencv/build/OpenCV-2.4.2/optimized/install/include

CFLAGS+=-Wimplicit -Werror=implicit-function-declaration -Werror=int-to-pointer-cast -I$(BASEDIR)/opencl -I$(BASEDIR)/core

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
export LDLIBS
export NO_INTEL_LIBS

all:
	+make -C core
	+make -C opencl
	+make -C base
	+make -C mlp
	+make -C OpenCV

clean:
	make -C core clean
	make -C opencl clean
	make -C base clean
	make -C mlp clean
	make -C OpenCV clean
# LuM end of file
