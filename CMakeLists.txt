# -*- sh -*-
# UjoImro, 2013

cmake_minimum_required (VERSION 2.8)

project(CarpBenchmark)

include(CheckIncludeFileCXX)
include(FindTBB.cmake)

if(NOT CMAKE_BUILD_TYPE)
    SET( CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build, options are: Release Debug" FORCE )
endif(NOT CMAKE_BUILD_TYPE)

string( COMPARE EQUAL ${CMAKE_HOST_SYSTEM_NAME} Linux THIS_IS_LINUX )

if (THIS_IS_LINUX)
    set( LIBEXT "so" )
    set( LIBPREFIX "lib" )
elseif (CMAKE_HOST_APPLE)
    set( LIBEXT "dylib" )
    set( LIBPREFIX "lib" )
elseif (CMAKE_HOST_WIN32)
    set( LIBEXT dll )
endif ()

set( POTENTIAL_LIBDIRS /usr/lib /usr/local/lib /opt/local/lib /usr/lib64 /usr/local/lib64 /opt/local/lib64 /Library )

################## PROCESSOR ####################
if(${CMAKE_SYSTEM_PROCESSOR} MATCHES amd64*)
   set(X86_64 1)
endif()

if(${CMAKE_SYSTEM_PROCESSOR} MATCHES x86_64*)
   set(X86_64 1)
endif()

if(${CMAKE_SYSTEM_PROCESSOR} MATCHES i686*)
   set(X86 1)
endif()

if(${CMAKE_SYSTEM_PROCESSOR} MATCHES i386*)
   set(X86 1)
endif()

if( (NOT X86_64) AND (${CMAKE_SYSTEM_PROCESSOR} MATCHES x86*) )
   set(X86 1)
endif()

if(CMAKE_CL_64 OR (${MSVC_C_ARCHITECTURE_ID} MATCHES x64))
   set(X86_64 1)
   unset(X86)
endif()

if(X86 AND NOT X86_64)
    message(STATUS "32-bit build detected.")
elseif(X86_64 AND NOT X86)
    message(STATUS "64-bit build detected.")
else()
    message(FATAL_ERROR  "Couldn't detect the correct bitness.")
endif()

if(FORCE_INTEL_COMPILER_SETTINGS OR CMAKE_C_COMPILER_ID STREQUAL "Intel")
    set(INTEL_COMPILER 1)
endif()

set(BUILD_SHARED_LIBS ON CACHE BOOL 
    "Build shared libraries (.dll/.so) instead of static ones (.lib/.a)")

if(X86)
    set(ENABLE_SSE OFF CACHE BOOL "Enable SSE/SSE2 instructions")
elseif(X86_64)
    set(ENABLE_SSE ON CACHE BOOL "Enable SSE/SSE2 instructions")
else()
    unset(ENABLE_SSE)
endif()

set(CMAKE_DEBUG_POSTFIX d)

if(X86)
    set( LIB_INSTALL_DIR lib )
elseif(X86_64)
    set( LIB_INSTALL_DIR lib64 )
endif()

################## XXD Compiler ####################
# The XXD Compiler compiles the OpenCL source code into C headers
find_program( XXD_COMPILER "xxd" )
if (NOT XXD_COMPILER)
    message( FATAL_ERROR "The xxd compiler not found. Please install the compiler before building.")
endif(NOT XXD_COMPILER)



######################### Includes ##########################
include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/include )
include_directories( ${CMAKE_BINARY_DIR} )

################### OPENCV ######################
# SET OpenCV PATHS
find_path(OpenCV_DIR "OpenCVConfig.cmake" DOC "Directory of the OpenCVConfig.cmake")
if(IS_DIRECTORY ${OpenCV_DIR})
    find_package(OpenCV REQUIRED)
    if (OpenCV_FOUND AND NOT (OpenCV_VERSION VERSION_LESS "2.4.5"))
        add_definitions(-D_CRT_SECURE_NO_WARNINGS)
    else()
        message(FATAL_ERROR "OpenCV version is less than expected. Found: ${OpenCV_VERSION} expected 2.4.5+")
    endif()

    set( CMAKE_REQUIRED_INCLUDES ${CMAKE_REQUIRED_INCLUDES} ${OpenCV_INCLUDE_DIRS} )

    CHECK_INCLUDE_FILE_CXX( opencv2/core/core.hpp HEADERTEST_OPENCV )
    if (NOT HEADERTEST_OPENCV)
	message( FATAL_ERROR "OpenCV Header does not Compile!")
    endif(NOT HEADERTEST_OPENCV)

else(IS_DIRECTORY ${OpenCV_DIR})
    message(FATAL_ERROR "OpenCV not found")
endif(IS_DIRECTORY ${OpenCV_DIR})

set( LIBS ${LIBS} ${OpenCV_LIBS} )

################### INTEL MKL ######################

SET( MKL_PATH "/opt/intel/mkl" CACHE PATH "The directory of the intel mkl library." )

option( WITH_MKL "Compiles the CLM module with the intel mkl library functions where appropriate" OFF )

if (WITH_MKL)
    if (LIBEXT STREQUAL "dll")
        set(MKL_LIB_POSTFIX "_dll.lib")
    else()
        set(MKL_LIB_POSTFIX ".${LIBEXT}")
    endif()

    if (NOT MKL_PATH)
	set( MKL_PATH /opt/intel/mkl CACHE PATH "The directory where then Intel MKL library is installed" )    
    endif( NOT MKL_PATH)

    if (NOT MKL_INCLUDE_PATH)
	set( MKL_INCLUDE_PATH ${MKL_PATH}/include CACHE PATH "The directory of the header files of the Intel MKL library." )
    endif()

    if (NOT MKL_LIB_PATH)
	if(X86_64)
     	    set( MKL_LIB_PATH "${MKL_PATH}/lib/intel64" CACHE PATH "The directory of the Intel MKL libraries." )
        elseif(X86)
            set( MKL_LIB_PATH "${MKL_PATH}/lib/ia32" CACHE PATH "The directory of the Intel MKL libraries." )
	endif()
    endif(NOT MKL_LIB_PATH)

    if (NOT CMAKE_HOST_WIN32)
        find_library( LIBM libm.${LIBEXT} PATHS ${POTENTIAL_LIBDIRS} ${MKL_LIB_PATH} )
	set( LIBS ${LIBS} ${LIBM} ) 
	
	find_library( LIB_PTHREAD libpthread.${LIBEXT} PATHS ${POTENTIAL_LIBDIRS} ${MKL_LIB_PATH} )
	set( LIBS ${LIBS} ${LIB_PTHREAD} )
    endif(NOT CMAKE_HOST_WIN32)

    if(X86_64)
        find_library( LIB_MKL_INTEL_LP64 ${LIBPREFIX}mkl_intel_lp64${MKL_LIB_POSTFIX} PATHS ${POTENTIAL_LIBDIRS} ${MKL_LIB_PATH} )
        set( LIBS ${LIBS} ${LIB_MKL_INTEL_LP64} )
    elseif(X86)
        find_library( LIB_MKL_INTEL_C ${LIBPREFIX}mkl_intel_c${MKL_LIB_POSTFIX} PATHS ${POTENTIAL_LIBDIRS} ${MKL_LIB_PATH} )
        set( LIBS ${LIBS} ${LIB_MKL_INTEL_C} )
    endif()

    find_library( LIB_SEQUENTIAL ${LIBPREFIX}mkl_sequential${MKL_LIB_POSTFIX} PATHS ${POTENTIAL_LIBDIRS} ${MKL_LIB_PATH} )
    set( LIBS ${LIBS} ${LIB_SEQUENTIAL} ) 

    find_library( LIB_MKL_CORE ${LIBPREFIX}mkl_core${MKL_LIB_POSTFIX} PATHS ${POTENTIAL_LIBDIRS} ${MKL_LIB_PATH} )
    set( LIBS ${LIBS} ${LIB_MKL_CORE} ) 

    link_directories( ${MKL_LIB_PATH} )    
    include_directories( ${MKL_INCLUDE_PATH} )

    add_definitions(-DWITH_MKL)

    ## No particular library or header file checking. It will fail during compilation, linking or first execution.
    ## If these tests are needed, we will set them up later.

    #include_directories( /home/ujoimro/opt/intel/composer_xe_2011_sp1.7.256/mkl/include )
    #link_directories( /home/ujoimro/opt/intel/composer_xe_2011_sp1.7.256/mkl/lib/intel64 )

    ## These are purely sequential
    #set( MKL_LIBS mkl_intel_lp64 mkl_sequential mkl_core pthread m )

    ## These are for openmp threading (it is generally of no use, as CLM already uses intel tbb)
    ## set( MKL_LIBS mkl_intel_lp64 mkl_intel_thread mkl_core iomp5 pthread m )
    ## this path was used for the intel openmp library 
    ## it hopefully will not be necessary
    ## link_directories( /home/ujoimro/opt/intel/composer_xe_2011_sp1.7.256/compiler/lib/intel64 )

endif(WITH_MKL)

################### TBB ######################
option( USE_TBB "Enables the Intel Thread Building Library" ON )

if (USE_TBB)
    add_definitions(-DWITH_TBB)

    if (NOT TBB_FOUND)
	message(FATAL_ERROR "TBB usage is set but TBB not found")
    endif(NOT TBB_FOUND)

    include_directories( TBB_INCLUDE_DIR )
    SET( LIBS ${LIBS} ${TBB_LIBRARIES} )
    # find_path(TBB_INCLUDE_PATH "tbb/tbb.h" PATHS "C:/TBB30/include" "/usr/include" )

    # if (THIS_IS_LINUX)
    #     if(X86_64)
    #         set( _TBB_LIB_PATH "${TBB_INCLUDE_PATH}/../lib64" )
    #     elseif(X86)
    #         set( _TBB_LIB_PATH "${TBB_INCLUDE_PATH}/../lib" )
    #     endif()

    # else (THIS_IS_LINUX)
    #     set(_TBB_LIB_PATH "${TBB_INCLUDE_PATH}/../lib")
    #     set(_TBB_BIN_PATH "${TBB_INCLUDE_PATH}/../bin")
    #     if(X86_64)
    #         set(_TBB_LIB_PATH "${_TBB_LIB_PATH}/intel64")
    #         set(_TBB_BIN_PATH "${_TBB_BIN_PATH}/intel64")
    #     elseif(X86)
    #         set(_TBB_LIB_PATH "${_TBB_LIB_PATH}/ia32")
    #         set(_TBB_BIN_PATH "${_TBB_BIN_PATH}/ia32")
    #     endif()
    # endif(THIS_IS_LINUX)

    # if(MSVC80)
    #     set(_TBB_LIB_PATH "${_TBB_LIB_PATH}/vc8")
    #     set(_TBB_BIN_PATH "${_TBB_BIN_PATH}/vc8")
    # elseif(MSVC90)
    #     set(_TBB_LIB_PATH "${_TBB_LIB_PATH}/vc9")
    #     set(_TBB_BIN_PATH "${_TBB_BIN_PATH}/vc9")
    # elseif(MSVC10)
    #     set(_TBB_LIB_PATH "${_TBB_LIB_PATH}/vc10")
    #     set(_TBB_BIN_PATH "${_TBB_BIN_PATH}/vc10")
    # elseif(MSVC11)
    #     set(_TBB_LIB_PATH "${_TBB_LIB_PATH}/vc11")
    #     set(_TBB_BIN_PATH "${_TBB_BIN_PATH}/vc11")
    # endif()

    # if (THIS_IS_LINUX)
    #     find_path(TBB_LIB_PATH libtbb.${LIBEXT} PATHS "${_TBB_LIB_PATH}")
    # else (THIS_IS_LINUX)
    #     find_path(TBB_LIB_PATH tbb.lib PATHS "${_TBB_LIB_PATH}")
    # endif (THIS_IS_LINUX)

    # if(IS_DIRECTORY ${TBB_INCLUDE_PATH} AND IS_DIRECTORY ${TBB_LIB_PATH})
    #     include_directories (${TBB_INCLUDE_PATH})
    #     link_directories ("${TBB_LIB_PATH}")

    #     if(NOT INTEL_COMPILER AND THIS_IS_LINUX)
    #         # assume TBB from repository
    #         set(TBB_LIBS
    #         debug tbb                 optimized tbb
    #         debug tbbmalloc           optimized tbbmalloc)
    #     else()
    #         # assume TBB from Intel compiler package
    #         set(TBB_LIBS
    #         debug tbb_debug           optimized tbb
    #         debug tbbmalloc_debug     optimized tbbmalloc)
    #     endif()

    # 	message( FATAL_ERROR ${TBB_LIBS} )
    #     set( LIBS ${LIBS} ${TBB_LIBS} )
    # else()
    #     message(FATAL_ERROR "TBB usage is set but TBB not found")
    # endif()

    # if (NOT WIN32)
    #     file(GLOB files "${TBB_LIB_PATH}/*tbb.so*")
    #     install( FILES ${files}
    #         DESTINATION ${LIB_INSTALL_DIR}
    #     )
    #     file(GLOB files "${TBB_LIB_PATH}/*tbbmalloc.so*")
    #     install( FILES ${files}
    #         DESTINATION ${LIB_INSTALL_DIR}
    #     )
    # endif(NOT WIN32)
endif(USE_TBB)


######################### BOOST ##########################
SET(BOOST_INCLUDEDIR "/usr/include" CACHE PATH "The directory of the boost headers. Remember in the case of /usr/local/include/boost, you specify /usr/local/include." )
if(CMAKE_HOST_WIN32)
    set(BOOST_LIBRARYDIR "/usr/lib" CACHE PATH "The directory of the boost libraries." )
endif()

set(BOOST_REQUIRED_COMPONENTS
        serialization
        iostreams
	thread
	system
)
set(BOOST_INCLUDE_CHECKS    
        boost/archive/xml_oarchive.hpp
        boost/archive/xml_iarchive.hpp
        boost/serialization/split_free.hpp
        boost/serialization/binary_object.hpp
        boost/serialization/serialization.hpp
)

if(BUILD_EXECUTABLES)
    set(BOOST_REQUIRED_COMPONENTS ${BOOST_REQUIRED_COMPONENTS}
            program_options
            chrono
            system
    )
    set(BOOST_INCLUDE_CHECKS ${BOOST_INCLUDE_CHECKS}
            boost/program_options.hpp
    )
endif()

if(WITH_TESTING)
    set(BOOST_REQUIRED_COMPONENTS ${BOOST_REQUIRED_COMPONENTS}
            program_options
            filesystem
            date_time
            system
    )
    set(BOOST_INCLUDE_CHECKS ${BOOST_INCLUDE_CHECKS}
            boost/program_options.hpp
    )

endif()

set(Boost_NO_SYSTEM_PATHS TRUE)
if(WIN32)
    set(Boost_USE_STATIC_RUNTIME OFF)
    set(Boost_USE_STATIC_LIBS ON)
endif()
if(FORCE_INTEL_COMPILER_SETTINGS AND WIN32)
    set(Boost_COMPILER "-iw")
    set(CMAKE_REQUIRED_FLAGS "/EHsc")
endif()

find_package(Boost COMPONENTS "${BOOST_REQUIRED_COMPONENTS}")

if(NOT Boost_FOUND)
    message( FATAL_ERROR "Boost libraries NOT FOUND!")
endif(NOT Boost_FOUND)

set( CMAKE_REQUIRED_INCLUDES ${CMAKE_REQUIRED_INCLUDES} ${Boost_INCLUDE_DIR} )
if(WIN32)
    set( CMAKE_REQUIRED_DEFINITIONS ${CMAKE_REQUIRED_DEFINITIONS} -DBOOST_ALL_NO_LIB )
endif()

CHECK_INCLUDE_FILE_CXX( opencv2/core/core.hpp HEADERTEST_OPENCV )
if (NOT HEADERTEST_OPENCV)
    message( FATAL_ERROR "OpenCV Header does not Compile!")
endif(NOT HEADERTEST_OPENCV)

foreach( header ${BOOST_INCLUDE_CHECKS} )
   string(REGEX REPLACE "[/\\.]" "_" header_name "${header}")
   CHECK_INCLUDE_FILE_CXX( ${header} HEADERTEST_${header_name} )
   if ( NOT HEADERTEST_${header_name} )
       message( FATAL_ERROR "Header " ${header} "not found!" )
   endif ( NOT HEADERTEST_${header_name} )
endforeach( header ${BOOST_INCLUDE_CHECKS} )

include_directories(${Boost_INCLUDE_DIR})
link_directories(${Boost_LIBRARY_DIRS})
set( LIBS ${LIBS} ${Boost_LIBRARIES})

if (THIS_IS_LINUX)
    file(GLOB files "${Boost_LIBRARY_DIRS}/libboost_*.so*")
    install( FILES ${files} DESTINATION ${LIB_INSTALL_DIR})
endif()


################## INTEL COMPILER SUPPORT LIBS #####################
if(INTEL_COMPILER)
    if(THIS_IS_LINUX)
        find_path(_INTEL_COMPILER_LIB_PATH "libimf.so" PATHS "${INTEL_COMPILER_LIB_PATH}" NO_DEFAULT_PATH)
        if (NOT _INTEL_COMPILER_LIB_PATH)
            message(FATAL_ERROR "Intel compiler selected but libimf.so not found.")
        endif()
      
        file(GLOB files "${_INTEL_COMPILER_LIB_PATH}/libimf.so*")
        install( FILES ${files} DESTINATION ${LIB_INSTALL_DIR})

        file(GLOB files "${_INTEL_COMPILER_LIB_PATH}/libsvml.so*")
        install( FILES ${files} DESTINATION ${LIB_INSTALL_DIR})

        file(GLOB files "${_INTEL_COMPILER_LIB_PATH}/libintlc.so*")
        install( FILES ${files} DESTINATION ${LIB_INSTALL_DIR})

        file(GLOB files "${_INTEL_COMPILER_LIB_PATH}/libirng.so*")
        install( FILES ${files} DESTINATION ${LIB_INSTALL_DIR})
    endif()
endif()

################## COMPILER SETTINGS #####################
add_definitions(-DNOMINMAX)

if (CMAKE_COMPILER_IS_GNUCXX)
    message(STATUS "Add GCC specific compiler flags")
    if(ENABLE_SSE AND X86)
        set(EXTRA_FLAGS "-mfpmath=sse -msse2")
    endif()
    set(CMAKE_C_FLAGS "-Wall -Wstrict-aliasing=2 ${EXTRA_FLAGS}")
    set(CMAKE_C_FLAGS_DEBUG "-O0 -ggdb3 -DDEBUG -D_DEBUG")
    set(CMAKE_C_FLAGS_RELEASE "-O3 -DNDEBUG -fomit-frame-pointer -ffast-math")
    set(CMAKE_CXX_FLAGS "-std=c++0x -Wall -Wstrict-aliasing=2 ${EXTRA_FLAGS}")
    set(CMAKE_CXX_FLAGS_DEBUG "-O0 -ggdb3 -DDEBUG -D_DEBUG")
    set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -fomit-frame-pointer -ffast-math")
endif()

if (MSVC)
    message(STATUS "Add MSVC specific compiler flags")
	set( WINDOWS_STORE_SUPPORT OFF CACHE BOOL "Compile for Windows Store support" )
	if( WINDOWS_STORE_SUPPORT )
		#TODO: update this when there is a new version of CMake
		message( "CMake does not support the global \"Compile for Windows Store\" switch yet, you have to turn that on manually.")
		add_definitions(-DUNICODE -D_UNICODE)
		set(EXTRA_FLAGS "/DWINAPI_FAMILY=WINAPI_FAMILY_APP /D_WIN32_WINNT=0x0602 /D_SCL_SECURE_NO_WARNINGS /Y-")
	endif()

    if(ENABLE_SSE)
        if(X86)
            set(EXTRA_FLAGS "/arch:SSE2 /D__SSE2__")
        elseif(X86_64)
            set(EXTRA_FLAGS "/D__SSE2__")
        endif()
    endif()
    set(CMAKE_C_FLAGS "/MP /W3 ${EXTRA_FLAGS}")
    set(CMAKE_C_FLAGS_DEBUG "/Od /Zi /DDEBUG /D_DEBUG /nologo /RTC1 /MDd /EHsc")
    set(CMAKE_C_FLAGS_RELEASE "/O2 /Ob2 /DNDEBUG /MD /Gy /EHsc /Zi /fp:fast")
    set(CMAKE_CXX_FLAGS "/MP /W3 ${EXTRA_FLAGS}")
    set(CMAKE_CXX_FLAGS_DEBUG "/Od /Zi /DDEBUG /D_DEBUG /nologo /RTC1 /MDd /EHsc")
    set(CMAKE_CXX_FLAGS_RELEASE "/O2 /Ob2 /DNDEBUG /MD /Gy /EHsc /Zi /fp:fast")
    set(CMAKE_EXE_LINKER_FLAGS "/debug")
    set(CMAKE_SHARED_LINKER_FLAGS "/debug")
	if( WINDOWS_STORE_SUPPORT )
		set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /APPCONTAINER")
		set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} /APPCONTAINER")
	endif()
endif()

if (INTEL_COMPILER)
    if (THIS_IS_LINUX)
        message(STATUS "Add Intel specific compiler flags (Linux)")

        if(ENABLE_SSE AND X86)
            set(EXTRA_FLAGS "-axSSE2")
        endif()
        set(CMAKE_C_FLAGS "-Wall ${EXTRA_FLAGS}")
        set(CMAKE_C_FLAGS_DEBUG "-O0 -ggdb3 -DDEBUG -D_DEBUG")
        set(CMAKE_C_FLAGS_RELEASE "-O3 -DNDEBUG -fomit-frame-pointer")
        set(CMAKE_CXX_FLAGS "-std=c++11 -Wall ${EXTRA_FLAGS}")
        set(CMAKE_CXX_FLAGS_DEBUG "-O0 -ggdb3 -DDEBUG -D_DEBUG")
        set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -fomit-frame-pointer")
    else()
        message(STATUS "Add Intel specific compiler flags (Windows)")

        set(CMAKE_C_FLAGS "/MP /W3 ${EXTRA_FLAGS}")
        set(CMAKE_C_FLAGS_DEBUG "/Od /Zi /DDEBUG /D_DEBUG /nologo /RTC1 /MDd /EHsc")
        set(CMAKE_C_FLAGS_RELEASE "/O2 /Ob2 /DNDEBUG /MD /Gy /EHsc /Zi /fp:fast")
        set(CMAKE_CXX_FLAGS "/Qstd=c++11 /MP /W3 ${EXTRA_FLAGS}")
        set(CMAKE_CXX_FLAGS_DEBUG "/Od /Zi /DDEBUG /D_DEBUG /nologo /RTC1 /MDd /EHsc")
        set(CMAKE_CXX_FLAGS_RELEASE "/O2 /Ob2 /DNDEBUG /MD /Gy /EHsc /Zi /fp:fast")
        set(CMAKE_EXE_LINKER_FLAGS "/debug")
        set(CMAKE_SHARED_LINKER_FLAGS "/debug")
    endif()
endif(INTEL_COMPILER)

################## OpenCL ######################

find_library( LIBOPENCL libOpenCL.${LIBEXT} PATHS ${POTENTIAL_LIBDIRS} )
set( LIBS ${LIBS} ${LIBOPENCL} )

################## Gathering Files ######################

INCLUDE(CheckTypeSize)

FILE( GLOB_RECURSE lib_src_c *.c )
FILE( GLOB_RECURSE lib_src_cpp *.cpp )
FILE( GLOB_RECURSE lib_src_cl *.cl )

## Collecting the Sources from the Directory
SET( SOURCES ${lib_src_c} ${lib_src_cpp} )

SET( TESTS "" )
SET( LIBSOURCES "" )
SET( CLS "" )
foreach( q ${SOURCES} )
    if ( q MATCHES "test_.*" )
	SET( TESTS ${TESTS} ${q} )
    else ( q MATCHES "test_.*" )
	SET( LIBSOURCES ${LIBSOURCES} ${q} )
    endif ( q MATCHES "test_.*" )
endforeach( q ${SOURCES} )

foreach( q ${lib_src_cl} )
    string( REGEX MATCH [^/]*.c[^.]*$ filename ${q} )
    string( LENGTH ${q} full_len )
    string( LENGTH ${filename} file_len )
    math( EXPR path_len  ${full_len}-${file_len} )
    string( SUBSTRING ${q} 0 ${path_len} path )
    string( REGEX REPLACE .c[^.]*$ "" executable ${filename} )

    message( STATUS "path = " ${path} )

    add_custom_command( OUTPUT ${filename}.clh
                        COMMAND "xxd" "-i" ${filename} ${CMAKE_BINARY_DIR}"/"${executable}.clh
                        WORKING_DIRECTORY ${path} )
    SET( CLS ${CLS} ${filename}.clh )
endforeach( q ${SOURCES} )

message( STATUS "TESTS = " ${TESTS} )
message( STATUS "LIBSOURCES = " ${LIBSOURCES} )

add_library( carpbench STATIC ${LIBSOURCES} ${CLS} )
target_link_libraries( carpbench ${LIBS} )

foreach( q ${TESTS} )
    string( REGEX MATCH [^/]*.c[^.]*$ filename ${q} )
    string( REGEX REPLACE .c[^.]*$ "" executable ${filename} )
    add_executable( ${executable} ${q} )
    add_dependencies( ${executable} carpbench )
    target_link_libraries( ${executable} carpbench ${LIBS} )
endforeach( q ${TESTS} )




# LuM end of file
