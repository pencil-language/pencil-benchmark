PPCG_COMPILER=/home/davidrobi/ppcg/ppcg
XXD_COMPILER=xxd

BOOST_INCLUDE_DIR=/usr/include/
BOOST_LIB_DIR=/usr/lib/x86_64-linux-gnu/
BOOST_LIBS=-lboost_filesystem -lboost_serialization -lboost_system

OPENCV_INCLUDE_DIR=/usr/local/include/
OPENCV_LIB_DIR=/usr/local/lib/
OPENCV_LIBS=-lopencv_core -lopencv_imgproc -lopencv_ocl -lopencv_highgui

TBB_INCLUDE_DIR=/usr/include/
TBB_LIB_DIR=/usr/lib/
TBB_LIBS=-ltbb -ltbbmalloc

OPENCL_PREFIX=/usr/
OPENCL_INCLUDE=$(OPENCL_PREFIX)/include/
OPENCL_LIB_DIR=$(OPENCL_PREFIX)/lib/
OPENCL_LIB=-lOpenCL

PPCG_OPTIONS=--no-shared-memory -D__PENCIL__ --target=opencl

# Optimization Flags

EXTRA_FLAGS=-O3 -DNDEBUG -march=native -fomit-frame-pointer -fPIC -ffast-math -Wall
# EXTRA_FLAGS=-O3 -DNDEBUG -march=native -fomit-frame-pointer -fPIC -ffast-math -Wall -DPRINT_OPENCL_PROFILING_KERNEL_EXEC_TIME
CFLAGS=$(EXTRA_FLAGS) -std=c99 -Iinclude -Ibuild -I$(OPENCL_INCLUDE)
CXXFLAGS=$(EXTRA_FLAGS) -std=c++0x -Iinclude -Ibuild -I$(OPENCL_INCLUDE) -I$(OPENCV_INCLUDE_DIR) -I$(TBB_INCLUDE_DIR) -I$(BOOST_INCLUDE_DIR)
LDFLAGS=-L$(OPENCL_LIB_DIR) $(OPENCL_LIB) -L$(OPENCV_LIB_DIR) $(OPENCV_LIBS) -L$(TBB_LIB_DIR) $(TBB_LIBS) -Lbuild -L$(BOOST_LIB_DIR) $(BOOST_LIBS) -Wl,-rpath=$$ORIGIN:$(OPENCV_LIB_DIR) -Wl,-z,origin

all: all_test all_ppcg_test mlp_data

clean: 
	-rm -f build/*.cl build/*.clh  build/*.c build/*.o  build/*.h build/*.so build/*.csv build/ppcg_test_* build/test_* build/temp_output_file build/temp_time* build/log

mlp_data: build/pool/response_dumps.xml

build/pool/response_dumps.xml:
	cd build/pool/; 7za e -y response_dumps.xml.7z


## Common Library
CL_SOURCES= ./gaussian/filter_sep_row.cl ./gaussian/filter_sep_col.cl ./cvt_color/cvt_color.cl ./filter2D/imgproc_convolve.cl ./dilate/filtering_morph.cl ./mlp/mlp_impl.cl ./mlp/operators.cl ./warpAffine/imgproc_warpAffine.cl ./resize/imgproc_resize.cl
LIB_SOURCES= 
PENCIL_SOURCES= ./gaussian/gaussian.pencil.c ./cvt_color/cvt_color.pencil.c ./filter2D/filter2D.pencil.c ./dilate/dilate.pencil.c ./warpAffine/warpAffine.pencil.c ./resize/resize.pencil.c ./hog/hog.pencil.c
MLP_SOURCES=./mlp/serialization.cpp ./mlp/allocator.cpp

## OpenCL Sources
all_opencl: build/filter_sep_row.clh build/filter_sep_col.clh build/cvt_color.clh build/imgproc_convolve.clh build/filtering_morph.clh build/mlp_impl.clh build/operators.clh build/imgproc_warpAffine.clh build/imgproc_resize.clh

build/filter_sep_row.clh: ./gaussian/filter_sep_row.cl
	(cd gaussian && $(XXD_COMPILER) -i filter_sep_row.cl ../build/filter_sep_row.clh)

build/filter_sep_col.clh: ./gaussian/filter_sep_col.cl
	(cd gaussian && $(XXD_COMPILER) -i filter_sep_col.cl ../build/filter_sep_col.clh)

build/cvt_color.clh: ./cvt_color/cvt_color.cl
	(cd cvt_color && $(XXD_COMPILER) -i cvt_color.cl ../build/cvt_color.clh)

build/imgproc_convolve.clh: ./filter2D/imgproc_convolve.cl
	(cd filter2D && $(XXD_COMPILER) -i imgproc_convolve.cl ../build/imgproc_convolve.clh)

build/filtering_morph.clh: ./dilate/filtering_morph.cl
	(cd dilate && $(XXD_COMPILER) -i filtering_morph.cl ../build/filtering_morph.clh)

build/mlp_impl.clh: ./mlp/mlp_impl.cl
	(cd mlp && $(XXD_COMPILER) -i mlp_impl.cl ../build/mlp_impl.clh)

build/operators.clh: ./mlp/operators.cl
	(cd mlp && $(XXD_COMPILER) -i operators.cl ../build/operators.clh)

build/imgproc_warpAffine.clh: ./warpAffine/imgproc_warpAffine.cl
	(cd warpAffine && $(XXD_COMPILER) -i imgproc_warpAffine.cl ../build/imgproc_warpAffine.clh)

build/imgproc_resize.clh: ./resize/imgproc_resize.cl
	(cd resize && $(XXD_COMPILER) -i imgproc_resize.cl ../build/imgproc_resize.clh)

## PENCIL GCC LIBRARY
build/libcarp_pencil.so: all_gcc_pencil_o
	$(CXX) -shared -o build/libcarp_pencil.so build/ocl_utilities.o  build/gaussian.pencil_as_c.o build/cvt_color.pencil_as_c.o build/filter2D.pencil_as_c.o build/dilate.pencil_as_c.o build/warpAffine.pencil_as_c.o build/resize.pencil_as_c.o build/mlp_impl.pencil_as_c.o build/hog.pencil_as_c.o $(LDFLAGS)

all_gcc_pencil_o: build/ocl_utilities.o build/gaussian.pencil_as_c.o build/cvt_color.pencil_as_c.o build/filter2D.pencil_as_c.o build/dilate.pencil_as_c.o build/warpAffine.pencil_as_c.o build/resize.pencil_as_c.o build/mlp_impl.pencil_as_c.o build/hog.pencil_as_c.o

build/gaussian.pencil_as_c.o: ./gaussian/gaussian.pencil.c
	$(CXX) -x c -c $(CFLAGS) ./gaussian/gaussian.pencil.c -o build/gaussian.pencil_as_c.o

build/cvt_color.pencil_as_c.o: ./cvt_color/cvt_color.pencil.c
	$(CXX) -x c -c $(CFLAGS) ./cvt_color/cvt_color.pencil.c -o build/cvt_color.pencil_as_c.o

build/filter2D.pencil_as_c.o: ./filter2D/filter2D.pencil.c
	$(CXX) -x c -c $(CFLAGS) ./filter2D/filter2D.pencil.c -o build/filter2D.pencil_as_c.o

build/dilate.pencil_as_c.o: ./dilate/dilate.pencil.c
	$(CXX) -x c -c $(CFLAGS) ./dilate/dilate.pencil.c -o build/dilate.pencil_as_c.o

build/warpAffine.pencil_as_c.o: ./warpAffine/warpAffine.pencil.c
	$(CXX) -x c -c $(CFLAGS) ./warpAffine/warpAffine.pencil.c -o build/warpAffine.pencil_as_c.o

build/resize.pencil_as_c.o: ./resize/resize.pencil.c
	$(CXX) -x c -c $(CFLAGS) ./resize/resize.pencil.c -o build/resize.pencil_as_c.o

build/mlp_impl.pencil_as_c.o: ./mlp/mlp_impl.pencil.c
	$(CXX) -x c -c $(CFLAGS) ./mlp/mlp_impl.pencil.c -o build/mlp_impl.pencil_as_c.o

build/hog.pencil_as_c.o: ./hog/hog.pencil.c
	$(CXX) -x c -c $(CFLAGS) ./hog/hog.pencil.c -o build/hog.pencil_as_c.o

build/ocl_utilities.o: ./base/ocl_utilities.c
	$(CXX) -x c -c $(CFLAGS) ./base/ocl_utilities.c -o build/ocl_utilities.o

## Standard Tests
all_test: build/test_gaussian build/test_cvt_color build/test_filter2D build/test_dilate build/test_mlp build/test_opencl_mlp build/test_gel_mlp build/test_warpAffine build/test_resize build/test_hog

build/test_gaussian: all_opencl ./gaussian/test_gaussian.cpp build/libcarp_pencil.so
	$(CXX) $(CXXFLAGS) ./gaussian/test_gaussian.cpp -o build/test_gaussian $(LDFLAGS) -lcarp_pencil 

build/test_cvt_color: all_opencl build/libcarp_pencil.so ./cvt_color/test_cvt_color.cpp
	$(CXX) $(CXXFLAGS) ./cvt_color/test_cvt_color.cpp -o build/test_cvt_color $(LDFLAGS) -lcarp_pencil 

build/test_filter2D: all_opencl build/libcarp_pencil.so ./filter2D/test_filter2D.cpp
	$(CXX) $(CXXFLAGS) ./filter2D/test_filter2D.cpp -o build/test_filter2D $(LDFLAGS) -lcarp_pencil 

build/test_dilate: all_opencl build/libcarp_pencil.so ./dilate/test_dilate.cpp
	$(CXX) $(CXXFLAGS) ./dilate/test_dilate.cpp -o build/test_dilate $(LDFLAGS) -lcarp_pencil 

build/test_mlp: all_opencl build/libcarp_pencil.so ./mlp/test_mlp.cpp $(MLP_SOURCES)
	$(CXX) $(CXXFLAGS) ./mlp/test_mlp.cpp $(MLP_SOURCES) -o build/test_mlp $(LDFLAGS) -lcarp_pencil

build/test_opencl_mlp: all_opencl build/libcarp_pencil.so ./mlp/test_opencl_mlp.cpp $(MLP_SOURCES)
	$(CXX) $(CXXFLAGS) ./mlp/test_opencl_mlp.cpp $(MLP_SOURCES) -o build/test_opencl_mlp $(LDFLAGS) -lcarp_pencil

build/test_gel_mlp: all_opencl build/libcarp_pencil.so ./mlp/test_gel_mlp.cpp $(MLP_SOURCES)
	$(CXX) $(CXXFLAGS) ./mlp/test_gel_mlp.cpp $(MLP_SOURCES) -o build/test_gel_mlp $(LDFLAGS) -lcarp_pencil 

build/test_warpAffine: all_opencl build/libcarp_pencil.so ./warpAffine/test_warpAffine.cpp 
	$(CXX) $(CXXFLAGS) ./warpAffine/test_warpAffine.cpp -o build/test_warpAffine $(LDFLAGS) -lcarp_pencil

build/test_resize: all_opencl build/libcarp_pencil.so ./resize/test_resize.cpp 
	$(CXX) $(CXXFLAGS) ./resize/test_resize.cpp -o build/test_resize $(LDFLAGS) -lcarp_pencil

build/test_hog: all_opencl build/libcarp_pencil.so ./hog/test_hog.cpp 
	$(CXX) $(CXXFLAGS) ./hog/test_hog.cpp -o build/test_hog $(LDFLAGS) -lcarp_pencil

## PPCG Compiled Source Files
all_pencil_source: build/gaussian.pencil_ppcg.c build/cvt_color.pencil_ppcg.c build/filter2D.pencil_ppcg.c build/dilate.pencil_ppcg.c build/warpAffine.pencil_ppcg.c build/resize.pencil_ppcg.c build/mlp_impl.pencil_ppcg.c  build/hog.pencil_ppcg.c

build/gaussian.pencil_ppcg.c: ./gaussian/gaussian.pencil.c
	$(PPCG_COMPILER) $(PPCG_OPTIONS) -o build/gaussian.pencil_ppcg.c gaussian/gaussian.pencil.c

build/cvt_color.pencil_ppcg.c: ./cvt_color/cvt_color.pencil.c
	$(PPCG_COMPILER) $(PPCG_OPTIONS) -o build/cvt_color.pencil_ppcg.c cvt_color/cvt_color.pencil.c

build/filter2D.pencil_ppcg.c: ./filter2D/filter2D.pencil.c
	$(PPCG_COMPILER) $(PPCG_OPTIONS) -o build/filter2D.pencil_ppcg.c filter2D/filter2D.pencil.c

build/dilate.pencil_ppcg.c: ./dilate/dilate.pencil.c
	$(PPCG_COMPILER) $(PPCG_OPTIONS) -o build/dilate.pencil_ppcg.c dilate/dilate.pencil.c

build/warpAffine.pencil_ppcg.c: ./warpAffine/warpAffine.pencil.c
	$(PPCG_COMPILER) $(PPCG_OPTIONS) -o build/warpAffine.pencil_ppcg.c warpAffine/warpAffine.pencil.c

build/resize.pencil_ppcg.c: ./resize/resize.pencil.c
	$(PPCG_COMPILER) $(PPCG_OPTIONS) -o build/resize.pencil_ppcg.c resize/resize.pencil.c

build/mlp_impl.pencil_ppcg.c: ./mlp/mlp_impl.pencil.c
	$(PPCG_COMPILER) $(PPCG_OPTIONS) -o build/mlp_impl.pencil_ppcg.c mlp/mlp_impl.pencil.c

build/hog.pencil_ppcg.c: ./hog/hog.pencil.c
	$(PPCG_COMPILER) $(PPCG_OPTIONS) -o build/hog.pencil_ppcg.c hog/hog.pencil.c

PPCG_INCLUDES=-I./gaussian -I./cvt_color -I./filter2D -I./dilate -I./warpAffine -I./resize -I./mlp -I./hog

build/warpAffine.pencil_ppcg.o: all_pencil_source
	$(CXX) -x c -c $(CFLAGS) $(PPCG_INCLUDES) build/warpAffine.pencil_ppcg.c -o build/warpAffine.pencil_ppcg.o

build/cvt_color.pencil_ppcg.o: all_pencil_source
	$(CXX) -x c -c $(CFLAGS) $(PPCG_INCLUDES) build/cvt_color.pencil_ppcg.c -o build/cvt_color.pencil_ppcg.o

build/dilate.pencil_ppcg.o: all_pencil_source
	$(CXX) -x c -c $(CFLAGS) $(PPCG_INCLUDES) build/dilate.pencil_ppcg.c -o build/dilate.pencil_ppcg.o

build/filter2D.pencil_ppcg.o: all_pencil_source
	$(CXX) -x c -c $(CFLAGS) $(PPCG_INCLUDES) build/filter2D.pencil_ppcg.c -o build/filter2D.pencil_ppcg.o

build/gaussian.pencil_ppcg.o: all_pencil_source
	$(CXX) -x c -c $(CFLAGS) $(PPCG_INCLUDES) build/gaussian.pencil_ppcg.c -o build/gaussian.pencil_ppcg.o

build/resize.pencil_ppcg.o: all_pencil_source
	$(CXX) -x c -c $(CFLAGS) $(PPCG_INCLUDES) build/resize.pencil_ppcg.c -o build/resize.pencil_ppcg.o

build/mlp_impl.pencil_ppcg.o: all_pencil_source
	$(CXX) -x c -c $(CFLAGS) $(PPCG_INCLUDES) build/mlp_impl.pencil_ppcg.c -o build/mlp_impl.pencil_ppcg.o

build/hog.pencil_ppcg.o: all_pencil_source
	$(CXX) -x c -c $(CFLAGS) $(PPCG_INCLUDES) build/hog.pencil_ppcg.c -o build/hog.pencil_ppcg.o

all_pencil_o: build/ocl_utilities.o build/warpAffine.pencil_ppcg.o build/cvt_color.pencil_ppcg.o build/dilate.pencil_ppcg.o build/filter2D.pencil_ppcg.o build/gaussian.pencil_ppcg.o build/resize.pencil_ppcg.o build/mlp_impl.pencil_ppcg.o build/hog.pencil_ppcg.o

build/libcarp_ppcg.so: all_pencil_o
	$(CXX) -shared -o  build/libcarp_ppcg.so build/ocl_utilities.o build/warpAffine.pencil_ppcg.o build/cvt_color.pencil_ppcg.o build/dilate.pencil_ppcg.o build/filter2D.pencil_ppcg.o build/gaussian.pencil_ppcg.o build/resize.pencil_ppcg.o build/mlp_impl.pencil_ppcg.o build/hog.pencil_ppcg.o $(LDFLAGS)

## PPCG Tests
all_ppcg_test: build/ppcg_test_gaussian build/ppcg_test_cvt_color build/ppcg_test_filter2D build/ppcg_test_dilate build/ppcg_test_mlp build/ppcg_test_opencl_mlp build/ppcg_test_gel_mlp build/ppcg_test_warpAffine build/ppcg_test_resize build/ppcg_test_hog

build/ppcg_test_gaussian: all_opencl build/libcarp_ppcg.so ./gaussian/test_gaussian.cpp
	$(CXX) $(CXXFLAGS) ./gaussian/test_gaussian.cpp -o build/ppcg_test_gaussian $(LDFLAGS) -lcarp_ppcg

build/ppcg_test_cvt_color: all_opencl build/libcarp_ppcg.so ./cvt_color/test_cvt_color.cpp
	$(CXX) $(CXXFLAGS) ./cvt_color/test_cvt_color.cpp -o build/ppcg_test_cvt_color $(LDFLAGS) -lcarp_ppcg

build/ppcg_test_filter2D: all_opencl build/libcarp_ppcg.so ./filter2D/test_filter2D.cpp
	$(CXX) $(CXXFLAGS) ./filter2D/test_filter2D.cpp -o build/ppcg_test_filter2D $(LDFLAGS) -lcarp_ppcg

build/ppcg_test_dilate: all_opencl build/libcarp_ppcg.so ./dilate/test_dilate.cpp
	$(CXX) $(CXXFLAGS) ./dilate/test_dilate.cpp -o build/ppcg_test_dilate $(LDFLAGS) -lcarp_ppcg

build/ppcg_test_mlp: all_opencl build/libcarp_ppcg.so ./mlp/test_mlp.cpp $(MLP_SOURCES)
	$(CXX) $(CXXFLAGS) ./mlp/test_mlp.cpp $(MLP_SOURCES) -o build/ppcg_test_mlp $(LDFLAGS) -lcarp_ppcg

build/ppcg_test_opencl_mlp: all_opencl build/libcarp_ppcg.so ./mlp/test_opencl_mlp.cpp $(MLP_SOURCES)
	$(CXX) $(CXXFLAGS) ./mlp/test_opencl_mlp.cpp $(MLP_SOURCES) -o build/ppcg_test_opencl_mlp $(LDFLAGS) -lcarp_ppcg

build/ppcg_test_gel_mlp: all_opencl build/libcarp_ppcg.so ./mlp/test_gel_mlp.cpp $(MLP_SOURCES)
	$(CXX) $(CXXFLAGS) ./mlp/test_gel_mlp.cpp $(MLP_SOURCES) -o build/ppcg_test_gel_mlp $(LDFLAGS) -lcarp_ppcg

build/ppcg_test_warpAffine: all_opencl build/libcarp_ppcg.so ./warpAffine/test_warpAffine.cpp 
	$(CXX) $(CXXFLAGS) ./warpAffine/test_warpAffine.cpp -o build/ppcg_test_warpAffine $(LDFLAGS) -lcarp_ppcg

build/ppcg_test_resize: all_opencl build/libcarp_ppcg.so ./resize/test_resize.cpp 
	$(CXX) $(CXXFLAGS) ./resize/test_resize.cpp -o build/ppcg_test_resize $(LDFLAGS) -lcarp_ppcg

build/ppcg_test_hog: all_opencl build/libcarp_ppcg.so ./hog/test_hog.cpp
	$(CXX) $(CXXFLAGS) ./hog/test_hog.cpp -o build/ppcg_test_hog $(LDFLAGS) -lcarp_ppcg
