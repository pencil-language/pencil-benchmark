cmake_minimum_required (VERSION 2.8.12)

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR})
set(CMAKE_DEBUG_POSTFIX d)
set(CMAKE_CONFIGURATION_TYPES Debug;Release)
include(PENCIL.cmake)

project(PENCIL_improc_benchmark)

######################### Compiler switches ##########################
add_compile_options(-march=native -flto -Wall -Wno-unknown-pragmas)

if (CMAKE_SYSTEM_PROCESSOR MATCHES "arm.*|ARM.*")
    add_compile_options(-mfpu=vfpv4)
else()
    add_compile_options(-mavx)
endif()
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c11")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")

########################## Dependencies ##########################
#directories: change in GUI if needed

find_package(OpenCV REQUIRED core imgproc highgui ocl)
find_package(TBB REQUIRED)

set(OPENCL_LIBRARY "/opt/AMDAPPSDK-2.9-1/lib/x86_64/libOpenCL.so.1" CACHE PATH "OpenCL runtime library")    #broken in FindOpenCL.cmake (upgrade to CMake 3.1,. or use -DOPENCL_LIBRARY=...) - Add libmali.so for Mali as well!
find_package(OpenCL REQUIRED)

set(Boost_USE_STATIC_LIBS       OFF)
set(Boost_USE_MULTITHREADED      ON)
set(Boost_USE_STATIC_RUNTIME    OFF)
set(Boost_ADDITIONAL_VERSIONS "1.54" "1.54.0" "1.55" "1.55.0")
find_package(Boost 1.54.0 REQUIRED filesystem serialization system)

######################### Add project files ##########################

set( COMMON_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/include
                         ${OpenCV_INCLUDE_DIRS}
                         ${Boost_INCLUDE_DIRS}
                         ${TBB_INCLUDE_DIRS}
                         ${OPENCL_INCLUDE_DIRS}
                         ${PENCIL_INCLUDE_DIRS}
                         )
set( COMMON_LINK_LIBRARIES ${OpenCV_LIBRARIES}
                           ${Boost_LIBRARIES}
                           ${TBB_LIBRARIES}
                           ${PENCIL_LIBRARIES}
                           ${OPENCL_LIBRARIES}
                           )

pencil_wrap(DEST cvt_color  DIM 2 FILES cvt_color/cvt_color.pencil.c)
pencil_wrap(DEST dilate     DIM 2 FILES dilate/dilate.pencil.c)
pencil_wrap(DEST filter2D   DIM 2 FILES filter2D/filter2D.pencil.c)
pencil_wrap(DEST gaussian   DIM 1 FILES gaussian/gaussian.pencil.c)
pencil_wrap(DEST histogram  DIM 1 FILES histogram/histogram.pencil.c)
pencil_wrap(DEST hog        DIM 1 FILES hog/hog.pencil.c)
pencil_wrap(DEST mlp        DIM 2 FILES mlp/mlp_impl.pencil.c)
pencil_wrap(DEST resize     DIM 2 FILES resize/resize.pencil.c)
pencil_wrap(DEST warpAffine DIM 2 FILES warpAffine/warpAffine.pencil.c)

set(cvt_color_SOURCES  cvt_color/test_cvt_color.cpp   cvt_color/cvt_color.pencil.h   )
set(dilate_SOURCES     dilate/test_dilate.cpp         dilate/dilate.pencil.h         )
set(filter2D_SOURCES   filter2D/test_filter2D.cpp     filter2D/filter2D.pencil.h     )
set(gaussian_SOURCES   gaussian/test_gaussian.cpp     gaussian/gaussian.pencil.h     )
set(histogram_SOURCES  histogram/test_histogram.cpp   histogram/histogram.pencil.h   )
set(hog_SOURCES hog/HogDescriptor.cpp
                hog/test_hog.cpp
                hog/hog.pencil.h
                hog/HogDescriptor.h
                )
set(mlp_SOURCES mlp/allocator.cpp
                mlp/serialization.cpp
                mlp/test_mlp.cpp
                mlp/cltypes.h
                mlp/mlp_impl_pencil.h
                mlp/allocator.hpp
                mlp/bench_mlp.hpp
                mlp/memory.hpp
                mlp/mlp.hpp
                mlp/serialization.hpp
                mlp/GEL/histo_impl.h
                mlp/GEL/histo.h
                )
set(resize_SOURCES     resize/test_resize.cpp         resize/resize.pencil.h         )
set(warpAffine_SOURCES warpAffine/test_warpAffine.cpp warpAffine/warpAffine.pencil.h )

add_executable(test_cvt_color  ${cvt_color_SOURCES}  ${cvt_color_GEN_SOURCES}  )
add_executable(test_dilate     ${dilate_SOURCES}     ${dilate_GEN_SOURCES}     )
add_executable(test_filter2D   ${filter2D_SOURCES}   ${filter2D_GEN_SOURCES}   )
add_executable(test_gaussian   ${gaussian_SOURCES}   ${gaussian_GEN_SOURCES}   )
add_executable(test_histogram  ${histogram_SOURCES}  ${histogram_GEN_SOURCES}  )
add_executable(test_hog        ${hog_SOURCES}        ${hog_GEN_SOURCES}        )
add_executable(test_mlp        ${mlp_SOURCES}        ${mlp_GEN_SOURCES}        )
add_executable(test_resize     ${resize_SOURCES}     ${resize_GEN_SOURCES}     )
add_executable(test_warpAffine ${warpAffine_SOURCES} ${warpAffine_GEN_SOURCES} )

target_include_directories( test_cvt_color  PRIVATE ${COMMON_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR}/cvt_color  ${cvt_color_GEN_INCLUDE_DIRS}  )
target_include_directories( test_dilate     PRIVATE ${COMMON_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR}/dilate     ${dilate_GEN_INCLUDE_DIRS}     )
target_include_directories( test_filter2D   PRIVATE ${COMMON_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR}/filter2D   ${filter2D_GEN_INCLUDE_DIRS}   )
target_include_directories( test_gaussian   PRIVATE ${COMMON_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR}/gaussian   ${gaussian_GEN_INCLUDE_DIRS}   )
target_include_directories( test_histogram  PRIVATE ${COMMON_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR}/histogram  ${histogram_GEN_INCLUDE_DIRS}  )
target_include_directories( test_hog        PRIVATE ${COMMON_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR}/hog        ${hog_GEN_INCLUDE_DIRS}        )
target_include_directories( test_mlp        PRIVATE ${COMMON_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR}/mlp        ${mlp_GEN_INCLUDE_DIRS}        )
target_include_directories( test_resize     PRIVATE ${COMMON_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR}/resize     ${resize_GEN_INCLUDE_DIRS}     )
target_include_directories( test_warpAffine PRIVATE ${COMMON_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR}/warpAffine ${warpAffine_GEN_INCLUDE_DIRS} )

target_link_libraries( test_cvt_color  ${COMMON_LINK_LIBRARIES} )
target_link_libraries( test_dilate     ${COMMON_LINK_LIBRARIES} )
target_link_libraries( test_filter2D   ${COMMON_LINK_LIBRARIES} )
target_link_libraries( test_gaussian   ${COMMON_LINK_LIBRARIES} )
target_link_libraries( test_histogram  ${COMMON_LINK_LIBRARIES} )
target_link_libraries( test_hog        ${COMMON_LINK_LIBRARIES} )
target_link_libraries( test_mlp        ${COMMON_LINK_LIBRARIES} )
target_link_libraries( test_resize     ${COMMON_LINK_LIBRARIES} )
target_link_libraries( test_warpAffine ${COMMON_LINK_LIBRARIES} )

add_custom_command( TARGET test_hog PRE_BUILD COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CMAKE_CURRENT_SOURCE_DIR}/hog/hog.opencl.cl ${CMAKE_CURRENT_BINARY_DIR}/kernels/hog.opencl.cl)
add_custom_command( TARGET test_mlp PRE_BUILD COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CMAKE_CURRENT_SOURCE_DIR}/mlp/mlp_impl.cl   ${CMAKE_CURRENT_BINARY_DIR}/kernels/mlp_impl.cl  )
