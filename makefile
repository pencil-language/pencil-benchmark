# UjoImro, 2013

# ## absolute directory for the subprojects
BASEDIR=$(shell pwd)

# Use icc by default, but allow the user to it
ifeq ($(USE_GCC), 1)
	CXX=g++
	CC=gcc
	LD=ld
else
	CXX?=icc
	CC?=icc
	LD?=xiar
endif

CXXFLAGS=-std=c++11 -I$(BASEDIR)/opencl -I$(BASEDIR)/core -gcc-name=gcc-4.6 -I/opt/intel/composer_xe_2013/ipp/include
CFLAGS=-lOpenCL
#LDFLAGS=-lirc -lOpenCL -lboost_serialization -L/home/ujoimro/Inst/opencv/build/OpenCV-2.4.2/optimized/install/lib -lopencv_core
#LDFLAGS=-lirc -lOpenCL -lboost_serialization -L/home/ujoimro/Inst/opencv/build/opencv-2.4.5/release/install/lib -lopencv_core -lopencv_ocl -lopencv_imgproc -lopencv_flann -lopencv_highgui -lopencv_features2d -lopencv_objdetect -lopencv_video -ltbb -lipps -lippi -lmkl_core -lmkl_vml_mc3 -lmkl_vml_avx -lmkl_vml_def -lmkl_vml_mc2 -lmkl_vml_p4n -lmkl_vml_mc -lmkl_vml_cmpt -lmkl_vml_avx2 -lmkl_intel_lp64 -lmkl_sequential
LDFLAGS=-L/home/ujoimro/Inst/opencv/build/opencv-2.4.5/release/install/lib

LDLIBS=-lirc -lOpenCL -lopencv_core -lopencv_ocl -lopencv_imgproc -lopencv_flann -lopencv_highgui -lopencv_features2d -lopencv_objdetect -lopencv_video -lboost_serialization -ltbb /usr/lib64/libm.so /usr/lib64/libpthread.so /opt/intel/mkl/lib/intel64/libmkl_intel_ilp64.so /opt/intel/mkl/lib/intel64/libmkl_sequential.so /opt/intel/mkl/lib/intel64/libmkl_core.so /opt/intel/ipp/lib/intel64/libippi.so /opt/intel/ipp/lib/intel64/libipps.so /opt/intel/ipp/lib/intel64/libippcc.so /opt/intel/ipp/lib/intel64/libippcore.so /opt/intel/ipp/lib/intel64/libippcv.so /opt/intel/composer_xe_2011_sp1.9.293/compiler/lib/intel64/libiomp5.so /opt/local/lib64/libboost_serialization.so /opt/local/lib64/libboost_iostreams.so /opt/local/lib64/libboost_thread.so /opt/local/lib64/libboost_system.so /opt/local/lib64/libboost_program_options.so /opt/local/lib64/libboost_chrono.so /opt/local/lib64/libboost_system.so /usr/lib64/libz.so 

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
export LDLIBS

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
