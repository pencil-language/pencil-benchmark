# UjoImro, 2013
# Experimental code for the CARP project

BUILD_DIR=build
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

PPCG_OPTIONS=--no-shared-memory -D__PENCIL__ 

# Optimization Flags

EXTRA_FLAGS=-O3 -DNDEBUG -march=native -fomit-frame-pointer -fPIC -ffast-math -Wall
# EXTRA_FLAGS=-O3 -DNDEBUG -march=native -fomit-frame-pointer -fPIC -ffast-math -Wall -DPRINT_OPENCL_PROFILING_KERNEL_EXEC_TIME
CFLAGS=$(EXTRA_FLAGS) -std=c99 -Iinclude -I$(BUILD_DIR) -I$(OPENCL_INCLUDE)
CXXFLAGS=$(EXTRA_FLAGS) -std=c++0x -Iinclude -I$(BUILD_DIR) -I$(OPENCL_INCLUDE) -I$(OPENCV_INCLUDE_DIR) -I$(TBB_INCLUDE_DIR) -I$(BOOST_INCLUDE_DIR)
LDFLAGS=-L$(OPENCL_LIB_DIR) $(OPENCL_LIB) -L$(OPENCV_LIB_DIR) $(OPENCV_LIBS) -L$(TBB_LIB_DIR) $(TBB_LIBS) -L$(BUILD_DIR) -L$(BOOST_LIB_DIR) $(BOOST_LIBS) -Wl,-rpath=$$ORIGIN:$(OPENCV_LIB_DIR) -Wl,-z,origin

all: all_test all_ppcg_test mlp_data

clean: 
	-rm -f $(BUILD_DIR)/*
	
mlp_data: $(BUILD_DIR)/pool/response_dumps.xml

$(BUILD_DIR)/pool/response_dumps.xml:
	cd $(BUILD_DIR)/pool/; 7za e -y response_dumps.xml.7z


## Common Library
CL_SOURCES= ./gaussian/filter_sep_row.cl ./gaussian/filter_sep_col.cl ./cvt_color/cvt_color.cl ./filter2D/imgproc_convolve.cl ./dilate/filtering_morph.cl ./mlp/mlp_impl.cl ./mlp/operators.cl ./warpAffine/imgproc_warpAffine.cl ./resize/imgproc_resize.cl
LIB_SOURCES= 
PENCIL_SOURCES= ./gaussian/gaussian.pencil.c ./cvt_color/cvt_color.pencil.c ./filter2D/filter2D.pencil.c ./dilate/dilate.pencil.c ./warpAffine/warpAffine.pencil.c ./resize/resize.pencil.c
MLP_SOURCES=./mlp/serialization.cpp ./mlp/allocator.cpp

## OpenCL Sources
all_opencl: $(BUILD_DIR)/filter_sep_row.clh $(BUILD_DIR)/filter_sep_col.clh $(BUILD_DIR)/cvt_color.clh $(BUILD_DIR)/imgproc_convolve.clh $(BUILD_DIR)/filtering_morph.clh $(BUILD_DIR)/mlp_impl.clh $(BUILD_DIR)/operators.clh $(BUILD_DIR)/imgproc_warpAffine.clh $(BUILD_DIR)/imgproc_resize.clh

$(BUILD_DIR)/filter_sep_row.clh: ./gaussian/filter_sep_row.cl
	ln -s ../gaussian/filter_sep_row.cl $(BUILD_DIR)/
	cd $(BUILD_DIR); $(XXD_COMPILER) -i filter_sep_row.cl filter_sep_row.clh

$(BUILD_DIR)/filter_sep_col.clh: ./gaussian/filter_sep_col.cl
	ln -s ../gaussian/filter_sep_col.cl $(BUILD_DIR)/
	cd $(BUILD_DIR); $(XXD_COMPILER) -i filter_sep_col.cl filter_sep_col.clh

$(BUILD_DIR)/cvt_color.clh: ./cvt_color/cvt_color.cl
	ln -s ../cvt_color/cvt_color.cl $(BUILD_DIR)/
	cd $(BUILD_DIR); $(XXD_COMPILER) -i cvt_color.cl cvt_color.clh

$(BUILD_DIR)/imgproc_convolve.clh: ./filter2D/imgproc_convolve.cl
	ln -s ../filter2D/imgproc_convolve.cl $(BUILD_DIR)/
	cd $(BUILD_DIR); $(XXD_COMPILER) -i imgproc_convolve.cl imgproc_convolve.clh

$(BUILD_DIR)/filtering_morph.clh: ./dilate/filtering_morph.cl
	ln -s ../dilate/filtering_morph.cl $(BUILD_DIR)/
	cd $(BUILD_DIR); $(XXD_COMPILER) -i filtering_morph.cl filtering_morph.clh

$(BUILD_DIR)/mlp_impl.clh: ./mlp/mlp_impl.cl
	ln -s ../mlp/mlp_impl.cl $(BUILD_DIR)/
	cd $(BUILD_DIR); $(XXD_COMPILER) -i mlp_impl.cl mlp_impl.clh

$(BUILD_DIR)/operators.clh: ./mlp/operators.cl
	ln -s ../mlp/operators.cl $(BUILD_DIR)/
	cd $(BUILD_DIR); $(XXD_COMPILER) -i operators.cl operators.clh

$(BUILD_DIR)/imgproc_warpAffine.clh: ./warpAffine/imgproc_warpAffine.cl
	ln -s ../warpAffine/imgproc_warpAffine.cl $(BUILD_DIR)/
	cd $(BUILD_DIR); $(XXD_COMPILER) -i imgproc_warpAffine.cl imgproc_warpAffine.clh

$(BUILD_DIR)/imgproc_resize.clh: ./resize/imgproc_resize.cl
	ln -s ../resize/imgproc_resize.cl $(BUILD_DIR)/
	cd $(BUILD_DIR); $(XXD_COMPILER) -i imgproc_resize.cl imgproc_resize.clh

## PENCIL GCC LIBRARY
$(BUILD_DIR)/libcarp_pencil.so: all_gcc_pencil_o
	$(CXX) -shared -o $(BUILD_DIR)/libcarp_pencil.so $(BUILD_DIR)/ocl_utilities.o  $(BUILD_DIR)/gaussian.pencil.o $(BUILD_DIR)/cvt_color.pencil.o $(BUILD_DIR)/filter2D.pencil.o $(BUILD_DIR)/dilate.pencil.o $(BUILD_DIR)/warpAffine.pencil.o $(BUILD_DIR)/resize.pencil.o $(BUILD_DIR)/mlp_impl.pencil.o $(LDFLAGS)

all_gcc_pencil_o: $(BUILD_DIR)/ocl_utilities.o $(BUILD_DIR)/gaussian.pencil.o $(BUILD_DIR)/cvt_color.pencil.o $(BUILD_DIR)/filter2D.pencil.o $(BUILD_DIR)/dilate.pencil.o $(BUILD_DIR)/warpAffine.pencil.o $(BUILD_DIR)/resize.pencil.o $(BUILD_DIR)/mlp_impl.pencil.o

$(BUILD_DIR)/gaussian.pencil.o: ./gaussian/gaussian.pencil.c
	$(CXX) -x c -c $(CFLAGS) ./gaussian/gaussian.pencil.c -o $(BUILD_DIR)/gaussian.pencil.o

$(BUILD_DIR)/cvt_color.pencil.o: ./cvt_color/cvt_color.pencil.c
	$(CXX) -x c -c $(CFLAGS) ./cvt_color/cvt_color.pencil.c -o $(BUILD_DIR)/cvt_color.pencil.o

$(BUILD_DIR)/filter2D.pencil.o: ./filter2D/filter2D.pencil.c
	$(CXX) -x c -c $(CFLAGS) ./filter2D/filter2D.pencil.c -o $(BUILD_DIR)/filter2D.pencil.o

$(BUILD_DIR)/dilate.pencil.o: ./dilate/dilate.pencil.c
	$(CXX) -x c -c $(CFLAGS) ./dilate/dilate.pencil.c -o $(BUILD_DIR)/dilate.pencil.o

$(BUILD_DIR)/warpAffine.pencil.o: ./warpAffine/warpAffine.pencil.c
	$(CXX) -x c -c $(CFLAGS) ./warpAffine/warpAffine.pencil.c -o $(BUILD_DIR)/warpAffine.pencil.o

$(BUILD_DIR)/resize.pencil.o: ./resize/resize.pencil.c
	$(CXX) -x c -c $(CFLAGS) ./resize/resize.pencil.c -o $(BUILD_DIR)/resize.pencil.o

$(BUILD_DIR)/mlp_impl.pencil.o: ./mlp/mlp_impl.pencil.c
	$(CXX) -x c -c $(CFLAGS) ./mlp/mlp_impl.pencil.c -o $(BUILD_DIR)/mlp_impl.pencil.o

$(BUILD_DIR)/ocl_utilities.o: ./base/ocl_utilities.c
	$(CXX) -x c -c $(CFLAGS) ./base/ocl_utilities.c -o $(BUILD_DIR)/ocl_utilities.o

## Standard Tests
all_test: $(BUILD_DIR)/test_gaussian $(BUILD_DIR)/test_cvt_color $(BUILD_DIR)/test_filter2D $(BUILD_DIR)/test_dilate $(BUILD_DIR)/test_mlp $(BUILD_DIR)/test_opencl_mlp $(BUILD_DIR)/test_gel_mlp $(BUILD_DIR)/test_warpAffine $(BUILD_DIR)/test_resize

$(BUILD_DIR)/test_gaussian: all_opencl ./gaussian/test_gaussian.cpp $(BUILD_DIR)/libcarp_pencil.so
	$(CXX) $(CXXFLAGS) ./gaussian/test_gaussian.cpp -o $(BUILD_DIR)/test_gaussian $(LDFLAGS) -lcarp_pencil 

$(BUILD_DIR)/test_cvt_color: all_opencl $(BUILD_DIR)/libcarp_pencil.so ./cvt_color/test_cvt_color.cpp
	$(CXX) $(CXXFLAGS) ./cvt_color/test_cvt_color.cpp -o $(BUILD_DIR)/test_cvt_color $(LDFLAGS) -lcarp_pencil 

$(BUILD_DIR)/test_filter2D: all_opencl $(BUILD_DIR)/libcarp_pencil.so ./filter2D/test_filter2D.cpp
	$(CXX) $(CXXFLAGS) ./filter2D/test_filter2D.cpp -o $(BUILD_DIR)/test_filter2D $(LDFLAGS) -lcarp_pencil 

$(BUILD_DIR)/test_dilate: all_opencl $(BUILD_DIR)/libcarp_pencil.so ./dilate/test_dilate.cpp
	$(CXX) $(CXXFLAGS) ./dilate/test_dilate.cpp -o $(BUILD_DIR)/test_dilate $(LDFLAGS) -lcarp_pencil 

$(BUILD_DIR)/test_mlp: all_opencl $(BUILD_DIR)/libcarp_pencil.so ./mlp/test_mlp.cpp $(MLP_SOURCES)
	$(CXX) $(CXXFLAGS) ./mlp/test_mlp.cpp $(MLP_SOURCES) -o $(BUILD_DIR)/test_mlp $(LDFLAGS) -lcarp_pencil

$(BUILD_DIR)/test_opencl_mlp: all_opencl $(BUILD_DIR)/libcarp_pencil.so ./mlp/test_opencl_mlp.cpp $(MLP_SOURCES)
	$(CXX) $(CXXFLAGS) ./mlp/test_opencl_mlp.cpp $(MLP_SOURCES) -o $(BUILD_DIR)/test_opencl_mlp $(LDFLAGS) -lcarp_pencil

$(BUILD_DIR)/test_gel_mlp: all_opencl $(BUILD_DIR)/libcarp_pencil.so ./mlp/test_gel_mlp.cpp $(MLP_SOURCES)
	$(CXX) $(CXXFLAGS) ./mlp/test_gel_mlp.cpp $(MLP_SOURCES) -o $(BUILD_DIR)/test_gel_mlp $(LDFLAGS) -lcarp_pencil 

$(BUILD_DIR)/test_warpAffine: all_opencl $(BUILD_DIR)/libcarp_pencil.so ./warpAffine/test_warpAffine.cpp 
	$(CXX) $(CXXFLAGS) ./warpAffine/test_warpAffine.cpp -o $(BUILD_DIR)/test_warpAffine $(LDFLAGS) -lcarp_pencil

$(BUILD_DIR)/test_resize: all_opencl $(BUILD_DIR)/libcarp_pencil.so ./resize/test_resize.cpp 
	$(CXX) $(CXXFLAGS) ./resize/test_resize.cpp -o $(BUILD_DIR)/test_resize $(LDFLAGS) -lcarp_pencil

## PPCG Compiled Source Files
all_pencil_source: $(BUILD_DIR)/gaussian.pencil_kernel.cl $(BUILD_DIR)/cvt_color.pencil_kernel.cl $(BUILD_DIR)/filter2D.pencil_kernel.cl $(BUILD_DIR)/dilate.pencil_kernel.cl $(BUILD_DIR)/warpAffine/warpAffine.pencil_kernel.cl $(BUILD_DIR)/resize/resize.pencil_kernel.cl $(BUILD_DIR)/mlp/mlp_impl.pencil_kernel.cl

$(BUILD_DIR)/gaussian.pencil_kernel.cl: ./gaussian/gaussian.pencil.c
	cd $(BUILD_DIR); $(PPCG_COMPILER) $(PPCG_OPTIONS) --target=opencl ../gaussian/gaussian.pencil.c

$(BUILD_DIR)/cvt_color.pencil_kernel.cl: ./cvt_color/cvt_color.pencil.c
	cd $(BUILD_DIR); $(PPCG_COMPILER) $(PPCG_OPTIONS) --target=opencl ../cvt_color/cvt_color.pencil.c

$(BUILD_DIR)/filter2D.pencil_kernel.cl: ./filter2D/filter2D.pencil.c
	cd $(BUILD_DIR); $(PPCG_COMPILER) $(PPCG_OPTIONS) --target=opencl ../filter2D/filter2D.pencil.c

$(BUILD_DIR)/dilate.pencil_kernel.cl: ./dilate/dilate.pencil.c
	cd $(BUILD_DIR); $(PPCG_COMPILER) $(PPCG_OPTIONS) --target=opencl ../dilate/dilate.pencil.c

$(BUILD_DIR)/warpAffine/warpAffine.pencil_kernel.cl: ./warpAffine/warpAffine.pencil.c
	cd $(BUILD_DIR); $(PPCG_COMPILER) $(PPCG_OPTIONS) --target=opencl ../warpAffine/warpAffine.pencil.c

$(BUILD_DIR)/resize/resize.pencil_kernel.cl: ./resize/resize.pencil.c
	cd $(BUILD_DIR); $(PPCG_COMPILER) $(PPCG_OPTIONS) --target=opencl ../resize/resize.pencil.c

$(BUILD_DIR)/mlp/mlp_impl.pencil_kernel.cl: ./mlp/mlp_impl.pencil.c
	cd $(BUILD_DIR); $(PPCG_COMPILER) $(PPCG_OPTIONS)--target=opencl ../mlp/mlp_impl.pencil.c

PPCG_INCLUDES=-I./gaussian -I./cvt_color -I./filter2D -I./dilate -I./warpAffine -I./resize -I./mlp

$(BUILD_DIR)/warpAffine.pencil_host.o: all_pencil_source
	$(CXX) -x c -c $(CFLAGS) $(PPCG_INCLUDES) $(BUILD_DIR)/warpAffine.pencil_host.c -o $(BUILD_DIR)/warpAffine.pencil_host.o

$(BUILD_DIR)/cvt_color.pencil_host.o: all_pencil_source
	$(CXX) -x c -c $(CFLAGS) $(PPCG_INCLUDES) $(BUILD_DIR)/cvt_color.pencil_host.c -o $(BUILD_DIR)/cvt_color.pencil_host.o

$(BUILD_DIR)/dilate.pencil_host.o: all_pencil_source
	$(CXX) -x c -c $(CFLAGS) $(PPCG_INCLUDES) $(BUILD_DIR)/dilate.pencil_host.c -o $(BUILD_DIR)/dilate.pencil_host.o

$(BUILD_DIR)/filter2D.pencil_host.o: all_pencil_source
	$(CXX) -x c -c $(CFLAGS) $(PPCG_INCLUDES) $(BUILD_DIR)/filter2D.pencil_host.c -o $(BUILD_DIR)/filter2D.pencil_host.o

$(BUILD_DIR)/gaussian.pencil_host.o: all_pencil_source
	$(CXX) -x c -c $(CFLAGS) $(PPCG_INCLUDES) $(BUILD_DIR)/gaussian.pencil_host.c -o $(BUILD_DIR)/gaussian.pencil_host.o

$(BUILD_DIR)/resize.pencil_host.o: all_pencil_source
	$(CXX) -x c -c $(CFLAGS) $(PPCG_INCLUDES) $(BUILD_DIR)/resize.pencil_host.c -o $(BUILD_DIR)/resize.pencil_host.o

$(BUILD_DIR)/mlp_impl.pencil_host.o: all_pencil_source
	$(CXX) -x c -c $(CFLAGS) $(PPCG_INCLUDES) $(BUILD_DIR)/mlp_impl.pencil_host.c -o $(BUILD_DIR)/mlp_impl.pencil_host.o

all_pencil_o: $(BUILD_DIR)/ocl_utilities.o $(BUILD_DIR)/warpAffine.pencil_host.o $(BUILD_DIR)/cvt_color.pencil_host.o $(BUILD_DIR)/dilate.pencil_host.o $(BUILD_DIR)/filter2D.pencil_host.o $(BUILD_DIR)/gaussian.pencil_host.o $(BUILD_DIR)/resize.pencil_host.o $(BUILD_DIR)/mlp_impl.pencil_host.o

$(BUILD_DIR)/libcarp_ppcg.so: all_pencil_o
	$(CXX) -shared -o  $(BUILD_DIR)/libcarp_ppcg.so $(BUILD_DIR)/ocl_utilities.o $(BUILD_DIR)/warpAffine.pencil_host.o $(BUILD_DIR)/cvt_color.pencil_host.o $(BUILD_DIR)/dilate.pencil_host.o $(BUILD_DIR)/filter2D.pencil_host.o $(BUILD_DIR)/gaussian.pencil_host.o $(BUILD_DIR)/resize.pencil_host.o $(BUILD_DIR)/mlp_impl.pencil_host.o $(LDFLAGS)

## PPCG Tests
all_ppcg_test: $(BUILD_DIR)/ppcg_test_gaussian $(BUILD_DIR)/ppcg_test_cvt_color $(BUILD_DIR)/ppcg_test_filter2D $(BUILD_DIR)/ppcg_test_dilate $(BUILD_DIR)/ppcg_test_mlp $(BUILD_DIR)/ppcg_test_opencl_mlp $(BUILD_DIR)/ppcg_test_gel_mlp $(BUILD_DIR)/ppcg_test_warpAffine $(BUILD_DIR)/ppcg_test_resize

$(BUILD_DIR)/ppcg_test_gaussian: all_opencl $(BUILD_DIR)/libcarp_ppcg.so ./gaussian/test_gaussian.cpp
	$(CXX) $(CXXFLAGS) ./gaussian/test_gaussian.cpp -o $(BUILD_DIR)/ppcg_test_gaussian $(LDFLAGS) -lcarp_ppcg

$(BUILD_DIR)/ppcg_test_cvt_color: all_opencl $(BUILD_DIR)/libcarp_ppcg.so ./cvt_color/test_cvt_color.cpp
	$(CXX) $(CXXFLAGS) ./cvt_color/test_cvt_color.cpp -o $(BUILD_DIR)/ppcg_test_cvt_color $(LDFLAGS) -lcarp_ppcg

$(BUILD_DIR)/ppcg_test_filter2D: all_opencl $(BUILD_DIR)/libcarp_ppcg.so ./filter2D/test_filter2D.cpp
	$(CXX) $(CXXFLAGS) ./filter2D/test_filter2D.cpp -o $(BUILD_DIR)/ppcg_test_filter2D $(LDFLAGS) -lcarp_ppcg

$(BUILD_DIR)/ppcg_test_dilate: all_opencl $(BUILD_DIR)/libcarp_ppcg.so ./dilate/test_dilate.cpp
	$(CXX) $(CXXFLAGS) ./dilate/test_dilate.cpp -o $(BUILD_DIR)/ppcg_test_dilate $(LDFLAGS) -lcarp_ppcg

$(BUILD_DIR)/ppcg_test_mlp: all_opencl $(BUILD_DIR)/libcarp_ppcg.so ./mlp/test_mlp.cpp $(MLP_SOURCES)
	$(CXX) $(CXXFLAGS) ./mlp/test_mlp.cpp $(MLP_SOURCES) -o $(BUILD_DIR)/ppcg_test_mlp $(LDFLAGS) -lcarp_ppcg

$(BUILD_DIR)/ppcg_test_opencl_mlp: all_opencl $(BUILD_DIR)/libcarp_ppcg.so ./mlp/test_opencl_mlp.cpp $(MLP_SOURCES)
	$(CXX) $(CXXFLAGS) ./mlp/test_opencl_mlp.cpp $(MLP_SOURCES) -o $(BUILD_DIR)/ppcg_test_opencl_mlp $(LDFLAGS) -lcarp_ppcg

$(BUILD_DIR)/ppcg_test_gel_mlp: all_opencl $(BUILD_DIR)/libcarp_ppcg.so ./mlp/test_gel_mlp.cpp $(MLP_SOURCES)
	$(CXX) $(CXXFLAGS) ./mlp/test_gel_mlp.cpp $(MLP_SOURCES) -o $(BUILD_DIR)/ppcg_test_gel_mlp $(LDFLAGS) -lcarp_ppcg

$(BUILD_DIR)/ppcg_test_warpAffine: all_opencl $(BUILD_DIR)/libcarp_ppcg.so ./warpAffine/test_warpAffine.cpp 
	$(CXX) $(CXXFLAGS) ./warpAffine/test_warpAffine.cpp -o $(BUILD_DIR)/ppcg_test_warpAffine $(LDFLAGS) -lcarp_ppcg

$(BUILD_DIR)/ppcg_test_resize: all_opencl $(BUILD_DIR)/libcarp_ppcg.so ./resize/test_resize.cpp 
	$(CXX) $(CXXFLAGS) ./resize/test_resize.cpp -o $(BUILD_DIR)/ppcg_test_resize $(LDFLAGS) -lcarp_ppcg
