## parameters, paths
## Edit these according to paths on your system.
PPCG_COMPILER=ppcg
XXD_COMPILER=xxd

BOOST_INCLUDE_DIR=/usr/include/
BOOST_LIB_DIR=/usr/lib/x86_64-linux-gnu/
BOOST_LIBS=-lboost_filesystem -lboost_serialization -lboost_system

OPENCV_PREFIX=/usr/local/
OPENCV_INCLUDE_DIR=$(OPENCV_PREFIX)include/
OPENCV_LIB_DIR=$(OPENCV_PREFIX)lib/
OPENCV_LIBS=-lopencv_core -lopencv_imgproc -lopencv_ocl -lopencv_highgui

TBB_INCLUDE_DIR=/usr/include/
TBB_LIB_DIR=/usr/lib/
TBB_LIBS=-ltbb -ltbbmalloc

OPENCL_PREFIX=/opt/AMDAPPSDK-2.9-1/
OPENCL_INCLUDE=$(OPENCL_PREFIX)include/
OPENCL_LIB_DIR=$(OPENCL_PREFIX)lib/
OPENCL_LIB=-lOpenCL
#OPENCL_LIB=-lOpenCL -lmali

OCL_UTILITIES=/home/davidrobi/CARP/ppcg/ocl_utilities.c
PENCIL_HEADERS=/home/davidrobi/CARP/pencil/etc/pencil_headers

# Optimization Flags

PPCG_OPTIONS=--no-shared-memory -D__PENCIL__ --target=opencl --sizes="{kernel[i]->block[16,16]}" -I$(PENCIL_HEADERS)
EXTRA_FLAGS=-O3 -DNDEBUG -march=native -fomit-frame-pointer -fPIC -ffast-math -Wall -Wno-unknown-pragmas -I. -I$(PENCIL_HEADERS)
#-DPRINT_OPENCL_PROFILING_KERNEL_EXEC_TIME

CFLAGS=$(EXTRA_FLAGS) -std=c1x -Iinclude -Ibuild -I$(OPENCL_INCLUDE)
CXXFLAGS=$(EXTRA_FLAGS) -std=c++0x -Iinclude -Ibuild -I$(OPENCL_INCLUDE) -I$(OPENCV_INCLUDE_DIR) -I$(TBB_INCLUDE_DIR) -I$(BOOST_INCLUDE_DIR)
LDFLAGS=-L$(OPENCL_LIB_DIR) $(OPENCL_LIB) -L$(OPENCV_LIB_DIR) $(OPENCV_LIBS) -L$(TBB_LIB_DIR) $(TBB_LIBS) -Lbuild -L$(BOOST_LIB_DIR) $(BOOST_LIBS) -Wl,-rpath=$$ORIGIN:$(OPENCV_LIB_DIR) -Wl,-z,origin ${EXTRA_OPENCL_LIBRARY}

#uncomment to compile with clang
#CXX=clang
#CXXFLAGS=$(EXTRA_FLAGS) -std=c++0x -Iinclude -Ibuild -I$(OPENCL_INCLUDE) -I$(OPENCV_INCLUDE_DIR) -I$(TBB_INCLUDE_DIR) -I$(BOOST_INCLUDE_DIR) -lstdc++ -lm

all: all_test all_ppcg_test mlp_data

all_test: build/test_gaussian build/test_cvt_color build/test_filter2D build/test_dilate build/test_mlp build/test_opencl_mlp build/test_gel_mlp build/test_warpAffine build/test_resize build/test_hog build/test_histogram

all_ppcg_test: build/ppcg_test_gaussian build/ppcg_test_cvt_color build/ppcg_test_filter2D build/ppcg_test_dilate build/ppcg_test_mlp build/ppcg_test_opencl_mlp build/ppcg_test_gel_mlp build/ppcg_test_warpAffine build/ppcg_test_resize build/ppcg_test_hog build/ppcg_test_histogram

all_pencil_source: build/gaussian.pencil_ppcg.c build/cvt_color.pencil_ppcg.c build/filter2D.pencil_ppcg.c build/dilate.pencil_ppcg.c build/warpAffine.pencil_ppcg.c build/resize.pencil_ppcg.c build/mlp_impl.pencil_ppcg.c build/hog.pencil_ppcg.c build/histogram.pencil_ppcg.c

clean: 
	@-rm -f build/*.cl build/*.clh  build/*.c build/*.o  build/*.h build/*.so build/*.csv build/ppcg_test_* build/test_* build/temp_output_file build/temp_time* build/log

mlp_data: build/pool/response_dumps.xml

build/pool/response_dumps.xml:
	@cd build/pool/; 7za e -y response_dumps.xml.7z

## OpenCL Sources
build/mlp_impl.clh: mlp/mlp_impl.cl
	@(cd mlp && $(XXD_COMPILER) -i mlp_impl.cl ../build/mlp_impl.clh)

build/operators.clh: mlp/operators.cl
	@(cd mlp && $(XXD_COMPILER) -i operators.cl ../build/operators.clh)

build/hog.clh: hog/hog.opencl.cl
	@(cd hog && $(XXD_COMPILER) -i hog.opencl.cl ../build/hog.clh)

## PENCIL OCL utilites
build/ocl_utilities.o:
	@$(CXX) -x c -c $(CFLAGS) $(OCL_UTILITIES) -o build/ocl_utilities.o

## PENCIL-as-c compile
build/cvt_color.pencil_as_c.o: cvt_color/cvt_color.pencil.c
	@$(CXX) -x c -c $(CFLAGS) cvt_color/cvt_color.pencil.c -o build/cvt_color.pencil_as_c.o

build/dilate.pencil_as_c.o: dilate/dilate.pencil.c
	@$(CXX) -x c -c $(CFLAGS) dilate/dilate.pencil.c -o build/dilate.pencil_as_c.o

build/filter2D.pencil_as_c.o: filter2D/filter2D.pencil.c
	@$(CXX) -x c -c $(CFLAGS) filter2D/filter2D.pencil.c -o build/filter2D.pencil_as_c.o

build/gaussian.pencil_as_c.o: gaussian/gaussian.pencil.c
	@$(CXX) -x c -c $(CFLAGS) gaussian/gaussian.pencil.c -o build/gaussian.pencil_as_c.o

build/histogram.pencil_as_c.o: histogram/histogram.pencil.c
	@$(CXX) -x c -c $(CFLAGS) histogram/histogram.pencil.c -o build/histogram.pencil_as_c.o

build/hog.pencil_as_c.o: hog/hog.pencil.c
	@$(CXX) -x c -c $(CFLAGS) hog/hog.pencil.c -o build/hog.pencil_as_c.o

build/mlp_impl.pencil_as_c.o: mlp/mlp_impl.pencil.c
	@$(CXX) -x c -c $(CFLAGS) mlp/mlp_impl.pencil.c -o build/mlp_impl.pencil_as_c.o

build/resize.pencil_as_c.o: resize/resize.pencil.c
	@$(CXX) -x c -c $(CFLAGS) resize/resize.pencil.c -o build/resize.pencil_as_c.o

build/warpAffine.pencil_as_c.o: warpAffine/warpAffine.pencil.c
	@$(CXX) -x c -c $(CFLAGS) warpAffine/warpAffine.pencil.c -o build/warpAffine.pencil_as_c.o



## PENCIL-as-c tests
build/test_cvt_color: cvt_color/test_cvt_color.cpp build/cvt_color.pencil_as_c.o build/ocl_utilities.o
	@$(CXX) $(CXXFLAGS) -o build/test_cvt_color cvt_color/test_cvt_color.cpp build/cvt_color.pencil_as_c.o build/ocl_utilities.o $(LDFLAGS)

build/test_dilate: dilate/test_dilate.cpp build/dilate.pencil_as_c.o build/ocl_utilities.o
	@$(CXX) $(CXXFLAGS) -o build/test_dilate dilate/test_dilate.cpp build/dilate.pencil_as_c.o build/ocl_utilities.o $(LDFLAGS)

build/test_filter2D: filter2D/test_filter2D.cpp build/filter2D.pencil_as_c.o build/ocl_utilities.o
	@$(CXX) $(CXXFLAGS) -o build/test_filter2D filter2D/test_filter2D.cpp build/filter2D.pencil_as_c.o build/ocl_utilities.o $(LDFLAGS)

build/test_gaussian: gaussian/test_gaussian.cpp build/gaussian.pencil_as_c.o build/ocl_utilities.o
	@$(CXX) $(CXXFLAGS) -o build/test_gaussian gaussian/test_gaussian.cpp build/gaussian.pencil_as_c.o build/ocl_utilities.o $(LDFLAGS)

build/test_histogram: histogram/test_histogram.cpp build/histogram.pencil_as_c.o build/ocl_utilities.o
	@$(CXX) $(CXXFLAGS) -o build/test_histogram histogram/test_histogram.cpp build/histogram.pencil_as_c.o build/ocl_utilities.o $(LDFLAGS)

build/test_hog: hog/test_hog.cpp build/hog.pencil_as_c.o build/ocl_utilities.o build/hog.clh
	@$(CXX) $(CXXFLAGS) -o build/test_hog hog/test_hog.cpp build/hog.pencil_as_c.o build/ocl_utilities.o $(LDFLAGS)

build/test_mlp: build/mlp_impl.clh build/operators.clh mlp/test_mlp.cpp mlp/serialization.cpp mlp/allocator.cpp build/mlp_impl.pencil_as_c.o build/ocl_utilities.o
	@$(CXX) $(CXXFLAGS) -o build/test_mlp mlp/test_mlp.cpp mlp/serialization.cpp mlp/allocator.cpp build/mlp_impl.pencil_as_c.o build/ocl_utilities.o $(LDFLAGS)

build/test_opencl_mlp: build/mlp_impl.clh build/operators.clh mlp/test_opencl_mlp.cpp mlp/serialization.cpp mlp/allocator.cpp build/mlp_impl.pencil_as_c.o build/ocl_utilities.o
	@$(CXX) $(CXXFLAGS) -o build/test_opencl_mlp mlp/test_opencl_mlp.cpp mlp/serialization.cpp mlp/allocator.cpp build/mlp_impl.pencil_as_c.o build/ocl_utilities.o $(LDFLAGS)

build/test_gel_mlp: build/mlp_impl.clh build/operators.clh mlp/test_gel_mlp.cpp mlp/serialization.cpp mlp/allocator.cpp build/mlp_impl.pencil_as_c.o build/ocl_utilities.o
	@$(CXX) $(CXXFLAGS) -o build/test_gel_mlp mlp/test_gel_mlp.cpp mlp/serialization.cpp mlp/allocator.cpp build/mlp_impl.pencil_as_c.o build/ocl_utilities.o $(LDFLAGS)

build/test_resize: resize/test_resize.cpp build/resize.pencil_as_c.o build/ocl_utilities.o
	@$(CXX) $(CXXFLAGS) -o build/test_resize resize/test_resize.cpp build/resize.pencil_as_c.o build/ocl_utilities.o $(LDFLAGS)

build/test_warpAffine: warpAffine/test_warpAffine.cpp build/warpAffine.pencil_as_c.o build/ocl_utilities.o
	@$(CXX) $(CXXFLAGS) -o build/test_warpAffine warpAffine/test_warpAffine.cpp build/warpAffine.pencil_as_c.o build/ocl_utilities.o $(LDFLAGS)



## PPCG compiled source files
build/cvt_color.pencil_ppcg.c: cvt_color/cvt_color.pencil.c
	@(cd build && $(PPCG_COMPILER) $(PPCG_OPTIONS) -o cvt_color.pencil_ppcg.c ../cvt_color/cvt_color.pencil.c )

build/dilate.pencil_ppcg.c: dilate/dilate.pencil.c
	@(cd build && $(PPCG_COMPILER) $(PPCG_OPTIONS) -o dilate.pencil_ppcg.c ../dilate/dilate.pencil.c )

build/filter2D.pencil_ppcg.c: filter2D/filter2D.pencil.c
	@(cd build && $(PPCG_COMPILER) $(PPCG_OPTIONS) -o filter2D.pencil_ppcg.c ../filter2D/filter2D.pencil.c )

build/gaussian.pencil_ppcg.c: gaussian/gaussian.pencil.c
	@(cd build && $(PPCG_COMPILER) $(PPCG_OPTIONS) -o gaussian.pencil_ppcg.c ../gaussian/gaussian.pencil.c )

build/histogram.pencil_ppcg.c: histogram/histogram.pencil.c
	@(cd build && $(PPCG_COMPILER) $(PPCG_OPTIONS) -o histogram.pencil_ppcg.c ../histogram/histogram.pencil.c )

build/hog.pencil_ppcg.c: hog/hog.pencil.c
	@(cd build && $(PPCG_COMPILER) $(PPCG_OPTIONS) -o hog.pencil_ppcg.c ../hog/hog.pencil.c )

build/mlp_impl.pencil_ppcg.c: mlp/mlp_impl.pencil.c
	@(cd build && $(PPCG_COMPILER) $(PPCG_OPTIONS) -o mlp_impl.pencil_ppcg.c ../mlp/mlp_impl.pencil.c )

build/resize.pencil_ppcg.c: resize/resize.pencil.c
	@(cd build && $(PPCG_COMPILER) $(PPCG_OPTIONS) -o resize.pencil_ppcg.c ../resize/resize.pencil.c )

build/warpAffine.pencil_ppcg.c: warpAffine/warpAffine.pencil.c
	@(cd build && $(PPCG_COMPILER) $(PPCG_OPTIONS) -o warpAffine.pencil_ppcg.c ../warpAffine/warpAffine.pencil.c )



##PENCIL host codes
build/cvt_color.pencil_ppcg.o: build/cvt_color.pencil_ppcg.c
	@$(CXX) -x c -c $(CFLAGS) -Icvt_color build/cvt_color.pencil_ppcg.c -o build/cvt_color.pencil_ppcg.o

build/dilate.pencil_ppcg.o: build/dilate.pencil_ppcg.c
	@$(CXX) -x c -c $(CFLAGS) -Idilate build/dilate.pencil_ppcg.c -o build/dilate.pencil_ppcg.o

build/filter2D.pencil_ppcg.o: build/filter2D.pencil_ppcg.c
	@$(CXX) -x c -c $(CFLAGS) -Ifilter2D build/filter2D.pencil_ppcg.c -o build/filter2D.pencil_ppcg.o

build/gaussian.pencil_ppcg.o: build/gaussian.pencil_ppcg.c
	@$(CXX) -x c -c $(CFLAGS) -Igaussian build/gaussian.pencil_ppcg.c -o build/gaussian.pencil_ppcg.o

build/histogram.pencil_ppcg.o: build/histogram.pencil_ppcg.c
	@$(CXX) -x c -c $(CFLAGS) -Ihistogram build/histogram.pencil_ppcg.c -o build/histogram.pencil_ppcg.o

build/hog.pencil_ppcg.o: build/hog.pencil_ppcg.c
	@$(CXX) -x c -c $(CFLAGS) -Ihog build/hog.pencil_ppcg.c -o build/hog.pencil_ppcg.o

build/mlp_impl.pencil_ppcg.o: build/mlp_impl.pencil_ppcg.c
	@$(CXX) -x c -c $(CFLAGS) -Imlp build/mlp_impl.pencil_ppcg.c -o build/mlp_impl.pencil_ppcg.o

build/resize.pencil_ppcg.o: build/resize.pencil_ppcg.c
	@$(CXX) -x c -c $(CFLAGS) -Iresize build/resize.pencil_ppcg.c -o build/resize.pencil_ppcg.o

build/warpAffine.pencil_ppcg.o: build/warpAffine.pencil_ppcg.c
	@$(CXX) -x c -c $(CFLAGS) -IwarpAffine build/warpAffine.pencil_ppcg.c -o build/warpAffine.pencil_ppcg.o



## PPCG tests
build/ppcg_test_cvt_color: cvt_color/test_cvt_color.cpp build/cvt_color.pencil_ppcg.o build/ocl_utilities.o
	@$(CXX) $(CXXFLAGS) -o build/ppcg_test_cvt_color cvt_color/test_cvt_color.cpp build/cvt_color.pencil_ppcg.o build/ocl_utilities.o $(LDFLAGS)

build/ppcg_test_dilate: dilate/test_dilate.cpp build/dilate.pencil_ppcg.o build/ocl_utilities.o
	@$(CXX) $(CXXFLAGS) -o build/ppcg_test_dilate dilate/test_dilate.cpp build/dilate.pencil_ppcg.o build/ocl_utilities.o $(LDFLAGS)

build/ppcg_test_filter2D: filter2D/test_filter2D.cpp build/filter2D.pencil_ppcg.o build/ocl_utilities.o
	@$(CXX) $(CXXFLAGS) -o build/ppcg_test_filter2D filter2D/test_filter2D.cpp build/filter2D.pencil_ppcg.o build/ocl_utilities.o $(LDFLAGS)

build/ppcg_test_gaussian: gaussian/test_gaussian.cpp build/gaussian.pencil_ppcg.o build/ocl_utilities.o
	@$(CXX) $(CXXFLAGS) -o build/ppcg_test_gaussian gaussian/test_gaussian.cpp build/gaussian.pencil_ppcg.o build/ocl_utilities.o $(LDFLAGS)

build/ppcg_test_histogram: histogram/test_histogram.cpp build/histogram.pencil_ppcg.o build/ocl_utilities.o
	@$(CXX) $(CXXFLAGS) -o build/ppcg_test_histogram histogram/test_histogram.cpp build/histogram.pencil_ppcg.o build/ocl_utilities.o $(LDFLAGS)

build/ppcg_test_hog: build/hog.clh hog/test_hog.cpp build/hog.pencil_ppcg.o build/ocl_utilities.o
	@$(CXX) $(CXXFLAGS) -o build/ppcg_test_hog hog/test_hog.cpp build/hog.pencil_ppcg.o build/ocl_utilities.o $(LDFLAGS)

build/ppcg_test_mlp: build/mlp_impl.clh build/operators.clh mlp/test_mlp.cpp mlp/serialization.cpp mlp/allocator.cpp build/mlp_impl.pencil_ppcg.o build/ocl_utilities.o
	@$(CXX) $(CXXFLAGS) -o build/ppcg_test_mlp mlp/test_mlp.cpp mlp/serialization.cpp mlp/allocator.cpp build/mlp_impl.pencil_ppcg.o build/ocl_utilities.o $(LDFLAGS)

build/ppcg_test_opencl_mlp: build/mlp_impl.clh build/operators.clh mlp/test_opencl_mlp.cpp mlp/serialization.cpp mlp/allocator.cpp build/mlp_impl.pencil_ppcg.o build/ocl_utilities.o
	@$(CXX) $(CXXFLAGS) -o build/ppcg_test_opencl_mlp mlp/test_opencl_mlp.cpp mlp/serialization.cpp mlp/allocator.cpp build/mlp_impl.pencil_ppcg.o build/ocl_utilities.o $(LDFLAGS)

build/ppcg_test_gel_mlp: build/mlp_impl.clh build/operators.clh mlp/test_gel_mlp.cpp mlp/serialization.cpp mlp/allocator.cpp build/mlp_impl.pencil_ppcg.o build/ocl_utilities.o
	@$(CXX) $(CXXFLAGS) -o build/ppcg_test_gel_mlp mlp/test_gel_mlp.cpp mlp/serialization.cpp mlp/allocator.cpp build/mlp_impl.pencil_ppcg.o build/ocl_utilities.o $(LDFLAGS)

build/ppcg_test_resize: resize/test_resize.cpp build/resize.pencil_ppcg.o build/ocl_utilities.o
	@$(CXX) $(CXXFLAGS) -o build/ppcg_test_resize resize/test_resize.cpp build/resize.pencil_ppcg.o build/ocl_utilities.o $(LDFLAGS)

build/ppcg_test_warpAffine: warpAffine/test_warpAffine.cpp build/warpAffine.pencil_ppcg.o build/ocl_utilities.o
	@$(CXX) $(CXXFLAGS) -o build/ppcg_test_warpAffine warpAffine/test_warpAffine.cpp build/warpAffine.pencil_ppcg.o build/ocl_utilities.o $(LDFLAGS)
