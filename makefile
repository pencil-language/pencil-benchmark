# UjoImro, 2013
# Experimental code for the CARP project

BUILD_DIR=build
PPCG_COMPILER=/home/davidrobi/ppcg/ppcg
XXD_COMPILER=xxd

BOOST_INCLUDE_DIR=/usr/include/
BOOST_LIB_DIR=/usr/lib/x86_64-linux-gnu/
BOOST_LIBS=-lboost_date_time -lboost_filesystem -lboost_iostreams -lboost_program_options -lboost_serialization -lboost_system -lboost_chrono

OPENCV_INCLUDE_DIR=/usr/local/include/
OPENCV_LIB_DIR=/usr/local/lib/
OPENCV_LIBS=-lopencv_calib3d -lopencv_flann -lopencv_legacy -lopencv_ocl -lopencv_ts -lopencv_contrib -lopencv_gpu -lopencv_ml -lopencv_photo -lopencv_video -lopencv_core -lopencv_highgui -lopencv_nonfree -lopencv_stitching -lopencv_videostab -lopencv_features2d -lopencv_imgproc -lopencv_objdetect -lopencv_superres

TBB_INCLUDE_DIR=/usr/include/
TBB_LIB_DIR=/usr/lib/
TBB_LIBS=-ltbb -ltbbmalloc

OPENCL_INCLUDE=/usr/include/
OPENCL_LIB_DIR=/usr/lib/
OPENCL_LIB=-lOpenCL

# Optimization Flags

EXTRA_FLAGS=-O3 -DNDEBUG -march=native -fomit-frame-pointer -fPIC -ffast-math -Wall -Wstrict-aliasing=2
CFLAGS=$(EXTRA_FLAGS) -std=c99 -Iinclude -I$(BUILD_DIR) -I$(OPENCL_INCLUDE)
CXXFLAGS=$(EXTRA_FLAGS) -std=c++0x -Iinclude -I$(BUILD_DIR) -I$(OPENCL_INCLUDE) -I$(OPENCV_INCLUDE_DIR) -I$(TBB_INCLUDE_DIR)
LDFLAGS=-L$(OPENCL_LIB_DIR) $(OPENCL_LIB) -L$(OPENCV_LIB_DIR) $(OPENCV_LIBS) -L$(TBB_LIB_DIR) $(TBB_LIBS) -L$(BUILD_DIR) $(BOOST_LIBS)

all: all_test all_ppcg_test

clean: 
	-rm -f -r $(BUILD_DIR)/*

## Common Library
CL_SOURCES= ./GaussianBlur/filter_sep_row.cl ./GaussianBlur/filter_sep_col.cl ./cvtColor/cvt_color.cl ./filter2D/imgproc_convolve.cl ./dilate/filtering_morph.cl ./cvIntegral/imgproc_integral_sum.cl ./mlp/mlp_impl.cl ./mlp/operators.cl ./boxFilter/filtering_boxFilter.cl ./warpAffine/imgproc_warpAffine.cl ./resize/imgproc_resize.cl ./base/nesting.cl ./base/local.cl ./base/color.cl 
LIB_SOURCES= 
PENCIL_SOURCES= ./GaussianBlur/gaussian.pencil.c ./cvtColor/cvt_color.pencil.c ./filter2D/filter2D.pencil.c ./dilate/dilate.pencil.c ./warpAffine/affine.pencil.c ./resize/resize.pencil.c
MLP_SOURCES=./mlp/serialization.cpp ./mlp/allocator.cpp ./mlp/GEL/linalg.cpp

## OpenCL Sources
all_opencl: $(BUILD_DIR)/filter_sep_row.clh $(BUILD_DIR)/filter_sep_col.clh $(BUILD_DIR)/cvt_color.clh $(BUILD_DIR)/imgproc_convolve.clh $(BUILD_DIR)/filtering_morph.clh $(BUILD_DIR)/imgproc_integral_sum.clh $(BUILD_DIR)/mlp_impl.clh $(BUILD_DIR)/operators.clh $(BUILD_DIR)/filtering_boxFilter.clh $(BUILD_DIR)/imgproc_warpAffine.clh $(BUILD_DIR)/imgproc_resize.clh $(BUILD_DIR)/nesting.clh $(BUILD_DIR)/local.clh $(BUILD_DIR)/color.clh

$(BUILD_DIR)/filter_sep_row.clh: ./GaussianBlur/filter_sep_row.cl
	ln -s ../GaussianBlur/filter_sep_row.cl $(BUILD_DIR)/
	cd $(BUILD_DIR); $(XXD_COMPILER) -i filter_sep_row.cl filter_sep_row.clh

$(BUILD_DIR)/filter_sep_col.clh: ./GaussianBlur/filter_sep_col.cl
	ln -s ../GaussianBlur/filter_sep_col.cl $(BUILD_DIR)/
	cd $(BUILD_DIR); $(XXD_COMPILER) -i filter_sep_col.cl filter_sep_col.clh

$(BUILD_DIR)/cvt_color.clh: ./cvtColor/cvt_color.cl
	ln -s ../cvtColor/cvt_color.cl $(BUILD_DIR)/
	cd $(BUILD_DIR); $(XXD_COMPILER) -i cvt_color.cl cvt_color.clh

$(BUILD_DIR)/imgproc_convolve.clh: ./filter2D/imgproc_convolve.cl
	ln -s ../filter2D/imgproc_convolve.cl $(BUILD_DIR)/
	cd $(BUILD_DIR); $(XXD_COMPILER) -i imgproc_convolve.cl imgproc_convolve.clh

$(BUILD_DIR)/filtering_morph.clh: ./dilate/filtering_morph.cl
	ln -s ../dilate/filtering_morph.cl $(BUILD_DIR)/
	cd $(BUILD_DIR); $(XXD_COMPILER) -i filtering_morph.cl filtering_morph.clh

$(BUILD_DIR)/imgproc_integral_sum.clh: ./cvIntegral/imgproc_integral_sum.cl
	ln -s ../cvIntegral/imgproc_integral_sum.cl $(BUILD_DIR)/
	cd $(BUILD_DIR); $(XXD_COMPILER) -i imgproc_integral_sum.cl imgproc_integral_sum.clh

$(BUILD_DIR)/mlp_impl.clh: ./mlp/mlp_impl.cl
	ln -s ../mlp/mlp_impl.cl $(BUILD_DIR)/
	cd $(BUILD_DIR); $(XXD_COMPILER) -i mlp_impl.cl mlp_impl.clh

$(BUILD_DIR)/operators.clh: ./mlp/operators.cl
	ln -s ../mlp/operators.cl $(BUILD_DIR)/
	cd $(BUILD_DIR); $(XXD_COMPILER) -i operators.cl operators.clh

$(BUILD_DIR)/filtering_boxFilter.clh: ./boxFilter/filtering_boxFilter.cl
	ln -s ../boxFilter/filtering_boxFilter.cl $(BUILD_DIR)/
	cd $(BUILD_DIR); $(XXD_COMPILER) -i filtering_boxFilter.cl filtering_boxFilter.clh

$(BUILD_DIR)/imgproc_warpAffine.clh: ./warpAffine/imgproc_warpAffine.cl
	ln -s ../warpAffine/imgproc_warpAffine.cl $(BUILD_DIR)/
	cd $(BUILD_DIR); $(XXD_COMPILER) -i imgproc_warpAffine.cl imgproc_warpAffine.clh

$(BUILD_DIR)/imgproc_resize.clh: ./resize/imgproc_resize.cl
	ln -s ../resize/imgproc_resize.cl $(BUILD_DIR)/
	cd $(BUILD_DIR); $(XXD_COMPILER) -i imgproc_resize.cl imgproc_resize.clh

$(BUILD_DIR)/nesting.clh: ./base/nesting.cl
	ln -s ../base/nesting.cl $(BUILD_DIR)/
	cd $(BUILD_DIR); $(XXD_COMPILER) -i nesting.cl nesting.clh

$(BUILD_DIR)/local.clh: ./base/local.cl
	ln -s ../base/local.cl $(BUILD_DIR)/
	cd $(BUILD_DIR); $(XXD_COMPILER) -i local.cl local.clh

$(BUILD_DIR)/color.clh: ./base/color.cl
	ln -s ../base/color.cl $(BUILD_DIR)/
	cd $(BUILD_DIR); $(XXD_COMPILER) -i color.cl color.clh

## PENCIL GCC LIBRARY
$(BUILD_DIR)/libcarp_pencil.so: all_gcc_pencil_o
	$(CXX) -shared -o $(BUILD_DIR)/libcarp_pencil.so $(BUILD_DIR)/ocl_utilities.o  $(BUILD_DIR)/gaussian.pencil.o $(BUILD_DIR)/cvt_color.pencil.o $(BUILD_DIR)/filter2D.pencil.o $(BUILD_DIR)/dilate.pencil.o $(BUILD_DIR)/affine.pencil.o $(BUILD_DIR)/resize.pencil.o $(BUILD_DIR)/mlp_impl.pencil.o $(LDFLAGS)

all_gcc_pencil_o: $(BUILD_DIR)/ocl_utilities.o $(BUILD_DIR)/gaussian.pencil.o $(BUILD_DIR)/cvt_color.pencil.o $(BUILD_DIR)/filter2D.pencil.o $(BUILD_DIR)/dilate.pencil.o $(BUILD_DIR)/affine.pencil.o $(BUILD_DIR)/resize.pencil.o $(BUILD_DIR)/mlp_impl.pencil.o

$(BUILD_DIR)/gaussian.pencil.o: ./GaussianBlur/gaussian.pencil.c
	$(CXX) -x c -c $(CFLAGS) ./GaussianBlur/gaussian.pencil.c -o $(BUILD_DIR)/gaussian.pencil.o

$(BUILD_DIR)/cvt_color.pencil.o: ./cvtColor/cvt_color.pencil.c
	$(CXX) -x c -c $(CFLAGS) ./cvtColor/cvt_color.pencil.c -o $(BUILD_DIR)/cvt_color.pencil.o

$(BUILD_DIR)/filter2D.pencil.o: ./filter2D/filter2D.pencil.c
	$(CXX) -x c -c $(CFLAGS) ./filter2D/filter2D.pencil.c -o $(BUILD_DIR)/filter2D.pencil.o

$(BUILD_DIR)/dilate.pencil.o: ./dilate/dilate.pencil.c
	$(CXX) -x c -c $(CFLAGS) ./dilate/dilate.pencil.c -o $(BUILD_DIR)/dilate.pencil.o

$(BUILD_DIR)/affine.pencil.o: ./warpAffine/affine.pencil.c
	$(CXX) -x c -c $(CFLAGS) ./warpAffine/affine.pencil.c -o $(BUILD_DIR)/affine.pencil.o

$(BUILD_DIR)/resize.pencil.o: ./resize/resize.pencil.c
	$(CXX) -x c -c $(CFLAGS) ./resize/resize.pencil.c -o $(BUILD_DIR)/resize.pencil.o

$(BUILD_DIR)/mlp_impl.pencil.o: ./mlp/mlp_impl.pencil.c
	$(CXX) -x c -c $(CFLAGS) ./mlp/mlp_impl.pencil.c -o $(BUILD_DIR)/mlp_impl.pencil.o

$(BUILD_DIR)/ocl_utilities.o: ./base/ocl_utilities.c
	$(CXX) -x c -c $(CFLAGS) ./base/ocl_utilities.c -o $(BUILD_DIR)/ocl_utilities.o

## Standard Tests
all_test: $(BUILD_DIR)/test_gaussian $(BUILD_DIR)/test_cvtColor $(BUILD_DIR)/test_filter2D $(BUILD_DIR)/test_dilate $(BUILD_DIR)/test_integral $(BUILD_DIR)/test_mlp $(BUILD_DIR)/test_opencl_mlp $(BUILD_DIR)/test_gel_mlp $(BUILD_DIR)/test_OpenCV $(BUILD_DIR)/test_boxFilter $(BUILD_DIR)/test_affine $(BUILD_DIR)/test_resize $(BUILD_DIR)/test_allocation $(BUILD_DIR)/test_local $(BUILD_DIR)/test_nesting $(BUILD_DIR)/test_color

$(BUILD_DIR)/test_gaussian: all_opencl ./GaussianBlur/test_gaussian.cpp $(BUILD_DIR)/libcarp_pencil.so
	$(CXX) $(CXXFLAGS) ./GaussianBlur/test_gaussian.cpp -o $(BUILD_DIR)/test_gaussian $(LDFLAGS) -lcarp_pencil 

$(BUILD_DIR)/test_cvtColor: all_opencl $(BUILD_DIR)/libcarp_pencil.so ./cvtColor/test_cvtColor.cpp
	$(CXX) $(CXXFLAGS) ./cvtColor/test_cvtColor.cpp -o $(BUILD_DIR)/test_cvtColor $(LDFLAGS) -lcarp_pencil 

$(BUILD_DIR)/test_filter2D: all_opencl $(BUILD_DIR)/libcarp_pencil.so ./filter2D/test_filter2D.cpp
	$(CXX) $(CXXFLAGS) ./filter2D/test_filter2D.cpp -o $(BUILD_DIR)/test_filter2D $(LDFLAGS) -lcarp_pencil 

$(BUILD_DIR)/test_dilate: all_opencl $(BUILD_DIR)/libcarp_pencil.so ./dilate/test_dilate.cpp
	$(CXX) $(CXXFLAGS) ./dilate/test_dilate.cpp -o $(BUILD_DIR)/test_dilate $(LDFLAGS) -lcarp_pencil 

$(BUILD_DIR)/test_integral: all_opencl $(BUILD_DIR)/libcarp_pencil.so ./cvIntegral/test_integral.cpp 
	$(CXX) $(CXXFLAGS) ./cvIntegral/test_integral.cpp -o $(BUILD_DIR)/test_integral $(LDFLAGS) -lcarp_pencil 

$(BUILD_DIR)/test_mlp: all_opencl $(BUILD_DIR)/libcarp_pencil.so ./mlp/test_mlp.cpp $(MLP_SOURCES)
	$(CXX) $(CXXFLAGS) ./mlp/test_mlp.cpp $(MLP_SOURCES) -o $(BUILD_DIR)/test_mlp $(LDFLAGS) -lcarp_pencil

$(BUILD_DIR)/test_opencl_mlp: all_opencl $(BUILD_DIR)/libcarp_pencil.so ./mlp/test_opencl_mlp.cpp $(MLP_SOURCES)
	$(CXX) $(CXXFLAGS) ./mlp/test_opencl_mlp.cpp $(MLP_SOURCES) -o $(BUILD_DIR)/test_opencl_mlp $(LDFLAGS) -lcarp_pencil

$(BUILD_DIR)/test_gel_mlp: all_opencl $(BUILD_DIR)/libcarp_pencil.so ./mlp/test_gel_mlp.cpp $(MLP_SOURCES)
	$(CXX) $(CXXFLAGS) ./mlp/test_gel_mlp.cpp $(MLP_SOURCES) -o $(BUILD_DIR)/test_gel_mlp $(LDFLAGS) -lcarp_pencil 

$(BUILD_DIR)/test_OpenCV: all_opencl $(BUILD_DIR)/libcarp_pencil.so ./OpenCV/test_OpenCV.cpp 
	$(CXX) $(CXXFLAGS) ./OpenCV/test_OpenCV.cpp -o $(BUILD_DIR)/test_OpenCV $(LDFLAGS) -lcarp_pencil

$(BUILD_DIR)/test_boxFilter: all_opencl $(BUILD_DIR)/libcarp_pencil.so ./boxFilter/test_boxFilter.cpp 
	$(CXX) $(CXXFLAGS) ./boxFilter/test_boxFilter.cpp -o $(BUILD_DIR)/test_boxFilter $(LDFLAGS) -lcarp_pencil

$(BUILD_DIR)/test_affine: all_opencl $(BUILD_DIR)/libcarp_pencil.so ./warpAffine/test_affine.cpp 
	$(CXX) $(CXXFLAGS) ./warpAffine/test_affine.cpp -o $(BUILD_DIR)/test_affine $(LDFLAGS) -lcarp_pencil

$(BUILD_DIR)/test_resize: all_opencl $(BUILD_DIR)/libcarp_pencil.so ./resize/test_resize.cpp 
	$(CXX) $(CXXFLAGS) ./resize/test_resize.cpp -o $(BUILD_DIR)/test_resize $(LDFLAGS) -lcarp_pencil

$(BUILD_DIR)/test_allocation: all_opencl $(BUILD_DIR)/libcarp_pencil.so ./base/test_allocation.cpp 
	$(CXX) $(CXXFLAGS) ./base/test_allocation.cpp -o $(BUILD_DIR)/test_allocation $(LDFLAGS) -lcarp_pencil

$(BUILD_DIR)/test_local: all_opencl $(BUILD_DIR)/libcarp_pencil.so ./base/test_local.cpp 
	$(CXX) $(CXXFLAGS) ./base/test_local.cpp -o $(BUILD_DIR)/test_local $(LDFLAGS) -lcarp_pencil

$(BUILD_DIR)/test_nesting: all_opencl $(BUILD_DIR)/libcarp_pencil.so ./base/test_nesting.cpp 
	$(CXX) $(CXXFLAGS) ./base/test_nesting.cpp -o $(BUILD_DIR)/test_nesting $(LDFLAGS) -lcarp_pencil

$(BUILD_DIR)/test_color: all_opencl $(BUILD_DIR)/libcarp_pencil.so ./base/test_color.cpp 
	$(CXX) $(CXXFLAGS) ./base/test_color.cpp -o $(BUILD_DIR)/test_color $(LDFLAGS) -lcarp_pencil

## PPCG Compiled Source Files
all_pencil_source: $(BUILD_DIR)/gaussian.pencil_kernel.cl $(BUILD_DIR)/cvt_color.pencil_kernel.cl $(BUILD_DIR)/filter2D.pencil_kernel.cl $(BUILD_DIR)/dilate.pencil_kernel.cl $(BUILD_DIR)/warpAffine/affine.pencil_kernel.cl $(BUILD_DIR)/resize/resize.pencil_kernel.cl $(BUILD_DIR)/mlp/mlp_impl.pencil_kernel.cl

$(BUILD_DIR)/gaussian.pencil_kernel.cl: ./GaussianBlur/gaussian.pencil.c
	cd $(BUILD_DIR); $(PPCG_COMPILER) --target=opencl ../GaussianBlur/gaussian.pencil.c

$(BUILD_DIR)/cvt_color.pencil_kernel.cl: ./cvtColor/cvt_color.pencil.c
	cd $(BUILD_DIR); $(PPCG_COMPILER) --target=opencl ../cvtColor/cvt_color.pencil.c

$(BUILD_DIR)/filter2D.pencil_kernel.cl: ./filter2D/filter2D.pencil.c
	cd $(BUILD_DIR); $(PPCG_COMPILER) --target=opencl ../filter2D/filter2D.pencil.c

$(BUILD_DIR)/dilate.pencil_kernel.cl: ./dilate/dilate.pencil.c
	cd $(BUILD_DIR); $(PPCG_COMPILER) --target=opencl ../dilate/dilate.pencil.c

$(BUILD_DIR)/warpAffine/affine.pencil_kernel.cl: ./warpAffine/affine.pencil.c
	cd $(BUILD_DIR); $(PPCG_COMPILER) --target=opencl ../warpAffine/affine.pencil.c

$(BUILD_DIR)/resize/resize.pencil_kernel.cl: ./resize/resize.pencil.c
	cd $(BUILD_DIR); $(PPCG_COMPILER) --target=opencl ../resize/resize.pencil.c

$(BUILD_DIR)/mlp/mlp_impl.pencil_kernel.cl: ./mlp/mlp_impl.pencil.c
	cd $(BUILD_DIR); $(PPCG_COMPILER) --target=opencl ../mlp/mlp_impl.pencil.c

PPCG_INCLUDES=-I./GaussianBlur -I./cvtColor -I./filter2D -I./dilate -I./warpAffine -I./resize -I./mlp

$(BUILD_DIR)/affine.pencil_host.o: all_pencil_source
	$(CXX) -x c -c $(CFLAGS) $(PPCG_INCLUDES) $(BUILD_DIR)/affine.pencil_host.c -o $(BUILD_DIR)/affine.pencil_host.o

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

all_pencil_o: $(BUILD_DIR)/ocl_utilities.o $(BUILD_DIR)/affine.pencil_host.o $(BUILD_DIR)/cvt_color.pencil_host.o $(BUILD_DIR)/dilate.pencil_host.o $(BUILD_DIR)/filter2D.pencil_host.o $(BUILD_DIR)/gaussian.pencil_host.o $(BUILD_DIR)/resize.pencil_host.o $(BUILD_DIR)/mlp_impl.pencil_host.o

$(BUILD_DIR)/libcarp_ppcg.so: all_pencil_o
	$(CXX) -shared -o  $(BUILD_DIR)/libcarp_ppcg.so $(BUILD_DIR)/ocl_utilities.o $(BUILD_DIR)/affine.pencil_host.o $(BUILD_DIR)/cvt_color.pencil_host.o $(BUILD_DIR)/dilate.pencil_host.o $(BUILD_DIR)/filter2D.pencil_host.o $(BUILD_DIR)/gaussian.pencil_host.o $(BUILD_DIR)/resize.pencil_host.o $(BUILD_DIR)/mlp_impl.pencil_host.o $(LDFLAGS)

## PPCG Tests
all_ppcg_test: $(BUILD_DIR)/ppcg_test_gaussian $(BUILD_DIR)/ppcg_test_cvtColor $(BUILD_DIR)/ppcg_test_filter2D $(BUILD_DIR)/ppcg_test_dilate $(BUILD_DIR)/ppcg_test_integral $(BUILD_DIR)/ppcg_test_mlp $(BUILD_DIR)/ppcg_test_opencl_mlp $(BUILD_DIR)/ppcg_test_gel_mlp $(BUILD_DIR)/ppcg_test_OpenCV $(BUILD_DIR)/ppcg_test_boxFilter $(BUILD_DIR)/ppcg_test_affine $(BUILD_DIR)/ppcg_test_resize $(BUILD_DIR)/ppcg_test_allocation $(BUILD_DIR)/ppcg_test_local $(BUILD_DIR)/ppcg_test_nesting $(BUILD_DIR)/ppcg_test_color

$(BUILD_DIR)/ppcg_test_gaussian: all_opencl $(BUILD_DIR)/libcarp_ppcg.so ./GaussianBlur/test_gaussian.cpp
	$(CXX) $(CXXFLAGS) ./GaussianBlur/test_gaussian.cpp -o $(BUILD_DIR)/ppcg_test_gaussian $(LDFLAGS) -lcarp_ppcg

$(BUILD_DIR)/ppcg_test_cvtColor: all_opencl $(BUILD_DIR)/libcarp_ppcg.so ./cvtColor/test_cvtColor.cpp
	$(CXX) $(CXXFLAGS) ./cvtColor/test_cvtColor.cpp -o $(BUILD_DIR)/ppcg_test_cvtColor $(LDFLAGS) -lcarp_ppcg

$(BUILD_DIR)/ppcg_test_filter2D: all_opencl $(BUILD_DIR)/libcarp_ppcg.so ./filter2D/test_filter2D.cpp
	$(CXX) $(CXXFLAGS) ./filter2D/test_filter2D.cpp -o $(BUILD_DIR)/ppcg_test_filter2D $(LDFLAGS) -lcarp_ppcg

$(BUILD_DIR)/ppcg_test_dilate: all_opencl $(BUILD_DIR)/libcarp_ppcg.so ./dilate/test_dilate.cpp
	$(CXX) $(CXXFLAGS) ./dilate/test_dilate.cpp -o $(BUILD_DIR)/ppcg_test_dilate $(LDFLAGS) -lcarp_ppcg

$(BUILD_DIR)/ppcg_test_integral: all_opencl $(BUILD_DIR)/libcarp_ppcg.so ./cvIntegral/test_integral.cpp 
	$(CXX) $(CXXFLAGS) ./cvIntegral/test_integral.cpp -o $(BUILD_DIR)/ppcg_test_integral $(LDFLAGS) -lcarp_ppcg

$(BUILD_DIR)/ppcg_test_mlp: all_opencl $(BUILD_DIR)/libcarp_ppcg.so ./mlp/test_mlp.cpp $(MLP_SOURCES)
	$(CXX) $(CXXFLAGS) ./mlp/test_mlp.cpp $(MLP_SOURCES) -o $(BUILD_DIR)/ppcg_test_mlp $(LDFLAGS) -lcarp_ppcg

$(BUILD_DIR)/ppcg_test_opencl_mlp: all_opencl $(BUILD_DIR)/libcarp_ppcg.so ./mlp/test_opencl_mlp.cpp $(MLP_SOURCES)
	$(CXX) $(CXXFLAGS) ./mlp/test_opencl_mlp.cpp $(MLP_SOURCES) -o $(BUILD_DIR)/ppcg_test_opencl_mlp $(LDFLAGS) -lcarp_ppcg

$(BUILD_DIR)/ppcg_test_gel_mlp: all_opencl $(BUILD_DIR)/libcarp_ppcg.so ./mlp/test_gel_mlp.cpp $(MLP_SOURCES)
	$(CXX) $(CXXFLAGS) ./mlp/test_gel_mlp.cpp $(MLP_SOURCES) -o $(BUILD_DIR)/ppcg_test_gel_mlp $(LDFLAGS) -lcarp_ppcg

$(BUILD_DIR)/ppcg_test_OpenCV: all_opencl $(BUILD_DIR)/libcarp_ppcg.so ./OpenCV/test_OpenCV.cpp 
	$(CXX) $(CXXFLAGS) ./OpenCV/test_OpenCV.cpp -o $(BUILD_DIR)/ppcg_test_OpenCV $(LDFLAGS) -lcarp_ppcg

$(BUILD_DIR)/ppcg_test_boxFilter: all_opencl $(BUILD_DIR)/libcarp_ppcg.so ./boxFilter/test_boxFilter.cpp 
	$(CXX) $(CXXFLAGS) ./boxFilter/test_boxFilter.cpp -o $(BUILD_DIR)/ppcg_test_boxFilter $(LDFLAGS) -lcarp_ppcg

$(BUILD_DIR)/ppcg_test_affine: all_opencl $(BUILD_DIR)/libcarp_ppcg.so ./warpAffine/test_affine.cpp 
	$(CXX) $(CXXFLAGS) ./warpAffine/test_affine.cpp -o $(BUILD_DIR)/ppcg_test_affine $(LDFLAGS) -lcarp_ppcg

$(BUILD_DIR)/ppcg_test_resize: all_opencl $(BUILD_DIR)/libcarp_ppcg.so ./resize/test_resize.cpp 
	$(CXX) $(CXXFLAGS) ./resize/test_resize.cpp -o $(BUILD_DIR)/ppcg_test_resize $(LDFLAGS) -lcarp_ppcg

$(BUILD_DIR)/ppcg_test_allocation: all_opencl $(BUILD_DIR)/libcarp_ppcg.so ./base/test_allocation.cpp 
	$(CXX) $(CXXFLAGS) ./base/test_allocation.cpp -o $(BUILD_DIR)/ppcg_test_allocation $(LDFLAGS) -lcarp_ppcg

$(BUILD_DIR)/ppcg_test_local: all_opencl $(BUILD_DIR)/libcarp_ppcg.so ./base/test_local.cpp 
	$(CXX) $(CXXFLAGS) ./base/test_local.cpp -o $(BUILD_DIR)/ppcg_test_local $(LDFLAGS) -lcarp_ppcg

$(BUILD_DIR)/ppcg_test_nesting: all_opencl $(BUILD_DIR)/libcarp_ppcg.so ./base/test_nesting.cpp 
	$(CXX) $(CXXFLAGS) ./base/test_nesting.cpp -o $(BUILD_DIR)/ppcg_test_nesting $(LDFLAGS) -lcarp_ppcg

$(BUILD_DIR)/ppcg_test_color: all_opencl $(BUILD_DIR)/libcarp_ppcg.so ./base/test_color.cpp 
	$(CXX) $(CXXFLAGS) ./base/test_color.cpp -o $(BUILD_DIR)/ppcg_test_color $(LDFLAGS) -lcarp_ppcg


# LuM end of file
