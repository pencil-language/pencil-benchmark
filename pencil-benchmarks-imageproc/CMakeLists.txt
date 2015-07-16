cmake_minimum_required (VERSION 2.8.7)

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR})
set(CMAKE_DEBUG_POSTFIX d)
set(CMAKE_CONFIGURATION_TYPES Debug;Release)
include(PENCIL.cmake)

project(PENCIL_improc_benchmark)

######################### Compiler switches ##########################
set(COMMON_COMPILE_FLAGS "-march=native -flto -Wall -Wno-unknown-pragmas")

if (CMAKE_SYSTEM_PROCESSOR MATCHES "arm.*|ARM.*")
    set(COMMON_COMPILE_FLAGS "${COMMON_COMPILE_FLAGS} -mfpu=vfpv4")
else()
    set(COMMON_COMPILE_FLAGS "${COMMON_COMPILE_FLAGS} -mavx")
endif()
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${COMMON_COMPILE_FLAGS} -std=c99")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${COMMON_COMPILE_FLAGS} -std=c++0x")

########################## Dependencies ##########################
#directories: change in GUI if needed

find_package(OpenCV REQUIRED core imgproc highgui ocl)
find_package(OpenCL REQUIRED)

set(Boost_USE_STATIC_LIBS       OFF)
set(Boost_USE_MULTITHREADED      ON)
set(Boost_USE_STATIC_RUNTIME    OFF)
find_package(Boost REQUIRED filesystem system)

######################### Optional dependencies ##########################
find_package(TBB)
if (TBB_LIBRARY)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DWITH_TBB")
endif()

######################### Add project files ##########################

set( COMMON_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/include
                         ${OpenCV_INCLUDE_DIRS}
                         ${Boost_INCLUDE_DIRS}
                         ${OPENCL_INCLUDE_DIRS}
                         ${PENCIL_INCLUDE_DIRS}
                         )
set( COMMON_LINK_LIBRARIES ${OpenCV_LIBRARIES}
                           ${Boost_LIBRARIES}
                           ${PENCIL_LIBRARIES}
                           ${OPENCL_LIBRARIES}
                           )

set(PENCIL_FLAGS_cvt_color  "" CACHE STRING "PENCIL compilation flags for cvt_color  - This is additional to PENCIL_REQUIRED_FLAGS. If not set, PENCIL_DEFAULT_* flags are used.")
set(PENCIL_FLAGS_dilate     "" CACHE STRING "PENCIL compilation flags for dilate     - This is additional to PENCIL_REQUIRED_FLAGS. If not set, PENCIL_DEFAULT_* flags are used.")
set(PENCIL_FLAGS_filter2D   "" CACHE STRING "PENCIL compilation flags for filter2D   - This is additional to PENCIL_REQUIRED_FLAGS. If not set, PENCIL_DEFAULT_* flags are used.")
set(PENCIL_FLAGS_gaussian   "" CACHE STRING "PENCIL compilation flags for gaussian   - This is additional to PENCIL_REQUIRED_FLAGS. If not set, PENCIL_DEFAULT_* flags are used.")
set(PENCIL_FLAGS_histogram  "" CACHE STRING "PENCIL compilation flags for histogram  - This is additional to PENCIL_REQUIRED_FLAGS. If not set, PENCIL_DEFAULT_* flags are used.")
set(PENCIL_FLAGS_hog        "" CACHE STRING "PENCIL compilation flags for hog        - This is additional to PENCIL_REQUIRED_FLAGS. If not set, PENCIL_DEFAULT_* flags are used.")
set(PENCIL_FLAGS_resize     "" CACHE STRING "PENCIL compilation flags for resize     - This is additional to PENCIL_REQUIRED_FLAGS. If not set, PENCIL_DEFAULT_* flags are used.")
set(PENCIL_FLAGS_warpAffine "" CACHE STRING "PENCIL compilation flags for warpAffine - This is additional to PENCIL_REQUIRED_FLAGS. If not set, PENCIL_DEFAULT_* flags are used.")

pencil_wrap(DEST cvt_color  FLAGS ${PENCIL_FLAGS_cvt_color}  FILES cvt_color/cvt_color.pencil.c)
pencil_wrap(DEST dilate     FLAGS ${PENCIL_FLAGS_dilate}     FILES dilate/dilate.pencil.c)
pencil_wrap(DEST filter2D   FLAGS ${PENCIL_FLAGS_filter2D}   FILES filter2D/filter2D.pencil.c)
pencil_wrap(DEST gaussian   FLAGS ${PENCIL_FLAGS_gaussian}   FILES gaussian/gaussian.pencil.c)
pencil_wrap(DEST histogram  FLAGS ${PENCIL_FLAGS_histogram}  FILES histogram/histogram.pencil.c)
pencil_wrap(DEST hog        FLAGS ${PENCIL_FLAGS_hog}        FILES hog/hog.pencil.c)
pencil_wrap(DEST resize     FLAGS ${PENCIL_FLAGS_resize}     FILES resize/resize.pencil.c)
pencil_wrap(DEST warpAffine FLAGS ${PENCIL_FLAGS_warpAffine} FILES warpAffine/warpAffine.pencil.c)

set(cvt_color_SOURCES  cvt_color/test_cvt_color.cpp   cvt_color/cvt_color.pencil.h   )
set(dilate_SOURCES     dilate/test_dilate.cpp         dilate/dilate.pencil.h         )
set(filter2D_SOURCES   filter2D/test_filter2D.cpp     filter2D/filter2D.pencil.h     )
set(gaussian_SOURCES   gaussian/test_gaussian.cpp     gaussian/gaussian.pencil.h     )
set(hog_SOURCES hog/HogDescriptor.cpp
                hog/test_hog.cpp
                hog/hog.pencil.h
                hog/HogDescriptor.h
                )
set(histogram_SOURCES  histogram/test_histogram.cpp   histogram/histogram.pencil.h   )
set(resize_SOURCES     resize/test_resize.cpp         resize/resize.pencil.h         )
set(warpAffine_SOURCES warpAffine/test_warpAffine.cpp warpAffine/warpAffine.pencil.h )

add_executable(test_cvt_color  ${cvt_color_SOURCES}  ${cvt_color_GEN_SOURCES}  )
add_executable(test_dilate     ${dilate_SOURCES}     ${dilate_GEN_SOURCES}     )
add_executable(test_filter2D   ${filter2D_SOURCES}   ${filter2D_GEN_SOURCES}   )
add_executable(test_gaussian   ${gaussian_SOURCES}   ${gaussian_GEN_SOURCES}   )
add_executable(test_histogram  ${histogram_SOURCES}  ${histogram_GEN_SOURCES}  )
add_executable(test_hog        ${hog_SOURCES}        ${hog_GEN_SOURCES}        )
add_executable(test_resize     ${resize_SOURCES}     ${resize_GEN_SOURCES}     )
add_executable(test_warpAffine ${warpAffine_SOURCES} ${warpAffine_GEN_SOURCES} )

if(${CMAKE_VERSION} VERSION_LESS 2.8.11)
    include_directories( ${COMMON_INCLUDE_DIRS}
                         ${CMAKE_CURRENT_SOURCE_DIR}/cvt_color  ${cvt_color_GEN_INCLUDE_DIRS}
                         ${CMAKE_CURRENT_SOURCE_DIR}/dilate     ${dilate_GEN_INCLUDE_DIRS}
                         ${CMAKE_CURRENT_SOURCE_DIR}/filter2D   ${filter2D_GEN_INCLUDE_DIRS}
                         ${CMAKE_CURRENT_SOURCE_DIR}/gaussian   ${gaussian_GEN_INCLUDE_DIRS}
                         ${CMAKE_CURRENT_SOURCE_DIR}/histogram  ${histogram_GEN_INCLUDE_DIRS}
                         ${CMAKE_CURRENT_SOURCE_DIR}/hog        ${hog_GEN_INCLUDE_DIRS}        ${TBB_INCLUDE_DIRS}
                         ${CMAKE_CURRENT_SOURCE_DIR}/resize     ${resize_GEN_INCLUDE_DIRS}
                         ${CMAKE_CURRENT_SOURCE_DIR}/warpAffine ${warpAffine_GEN_INCLUDE_DIRS}
                       )
else()
    target_include_directories( test_cvt_color  PRIVATE ${COMMON_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR}/cvt_color  ${cvt_color_GEN_INCLUDE_DIRS}  )
    target_include_directories( test_dilate     PRIVATE ${COMMON_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR}/dilate     ${dilate_GEN_INCLUDE_DIRS}     )
    target_include_directories( test_filter2D   PRIVATE ${COMMON_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR}/filter2D   ${filter2D_GEN_INCLUDE_DIRS}   )
    target_include_directories( test_gaussian   PRIVATE ${COMMON_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR}/gaussian   ${gaussian_GEN_INCLUDE_DIRS}   )
    target_include_directories( test_histogram  PRIVATE ${COMMON_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR}/histogram  ${histogram_GEN_INCLUDE_DIRS}  )
    target_include_directories( test_hog        PRIVATE ${COMMON_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR}/hog        ${hog_GEN_INCLUDE_DIRS}        ${TBB_INCLUDE_DIRS})
    target_include_directories( test_resize     PRIVATE ${COMMON_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR}/resize     ${resize_GEN_INCLUDE_DIRS}     )
    target_include_directories( test_warpAffine PRIVATE ${COMMON_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR}/warpAffine ${warpAffine_GEN_INCLUDE_DIRS} )
endif()

target_link_libraries( test_cvt_color  ${COMMON_LINK_LIBRARIES} )
target_link_libraries( test_dilate     ${COMMON_LINK_LIBRARIES} )
target_link_libraries( test_filter2D   ${COMMON_LINK_LIBRARIES} )
target_link_libraries( test_gaussian   ${COMMON_LINK_LIBRARIES} )
target_link_libraries( test_histogram  ${COMMON_LINK_LIBRARIES} )
target_link_libraries( test_hog        ${COMMON_LINK_LIBRARIES} ${TBB_LIBRARIES})
target_link_libraries( test_resize     ${COMMON_LINK_LIBRARIES} )
target_link_libraries( test_warpAffine ${COMMON_LINK_LIBRARIES} )

add_custom_command( TARGET test_hog PRE_BUILD COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CMAKE_CURRENT_SOURCE_DIR}/hog/hog.opencl.cl ${CMAKE_CURRENT_BINARY_DIR}/hog.opencl.cl)
