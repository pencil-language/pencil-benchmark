# UjoImro, 2013

# ## absolute directory for the subprojects
BASEDIR=$(shell pwd)

#CXX=g++-4.6
#CC=gcc-4.6
#LD=ld

CXX=icc
CC=icc
LD=xiar
CXXFLAGS=-std=c++11 -I$(BASEDIR)/opencl -I$(BASEDIR)/core -gcc-name=gcc-4.6
CFLAGS=-lOpenCL
#LDFLAGS=-lirc -lOpenCL -lboost_serialization -L/home/ujoimro/Inst/opencv/build/OpenCV-2.4.2/optimized/install/lib -lopencv_core
LDFLAGS=-lirc -lOpenCL -lboost_serialization -L/home/ujoimro/Inst/opencv/build/opencv-2.4.5/release/install/lib -lopencv_core -lopencv_ocl -lopencv_imgproc -lopencv_flann -lopencv_highgui -lopencv_features2d -lopencv_objdetect -lopencv_video

# optimization flags
#CXXFLAGS+=-O3 -mtune=native
CXXFLAGS+=-O3 -xHOST -ipo

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

all:
	/opt/intel/composer_xe_2013.4.183/bin/compilervars.sh intel64
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
