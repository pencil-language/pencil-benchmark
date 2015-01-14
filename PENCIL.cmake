set(PENCIL_COMPILER     "ppcg"                 CACHE STRING "PENCIL compiler")
set(PENCIL_INCLUDE_DIRS "/usr/local/include"   CACHE PATH   "PENCIL library headers directory")
set(PENCIL_LIBRARIES    "libocl_pencil_opt.so" CACHE FILE   "PENCIL runtime library")

set(PENCIL_FLAGS          "--sizes=\"{kernel[i]->block[16,16]}\""                                                                                  CACHE STRING "PENCIL compilation flags. The overall number of threads cannot be more than the OpenCL device allows (ARM Mali: 64-256 (depends on kernel), AMD: 256, Intel HD: 512, nVidia Fermi: 4096)")
set(PENCIL_REQUIRED_FLAGS "-D__PENCIL__;--target=opencl;--opencl-include-file=${PENCIL_INCLUDE_DIRS}/pencil_opencl.h;--isl-ast-always-print-block" CACHE STRING "Required PENCIL compilation flags")
mark_as_advanced(PENCIL_REQUIRED_FLAGS)

function(pencil_wrap)
  cmake_parse_arguments(COMPILE_PENCIL "" "DEST" "FLAGS;FILES" ${ARGN})
  set(COMPILE_PENCIL_DEST_INCLUDE_DIRS "${COMPILE_PENCIL_DEST}_GEN_INCLUDE_DIRS")
  set(COMPILE_PENCIL_DEST_SOURCES "${COMPILE_PENCIL_DEST}_GEN_SOURCES")
  if("${COMPILE_PENCIL_FLAGS}" STREQUAL "")
    message(STATUS "Flags not set for ${COMPILE_PENCIL_DEST}, using global PENCIL flags")
    set(COMPILE_PENCIL_FLAGS ${PENCIL_FLAGS})
  endif()
  foreach(pencil_file ${COMPILE_PENCIL_FILES})
    get_filename_component(PENCIL_FILE_NAME ${pencil_file} NAME_WE)
    get_filename_component(PENCIL_FILE_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/${pencil_file} DIRECTORY)
    list(APPEND ${COMPILE_PENCIL_DEST_INCLUDE_DIRS} ${PENCIL_FILE_DIRECTORY})
    add_custom_command(OUTPUT ${PENCIL_FILE_NAME}.ppcg.c ${PENCIL_FILE_NAME}.ppcg_kernel.cl
                       COMMAND ${PENCIL_COMPILER}
                       ARGS ${PENCIL_REQUIRED_FLAGS} ${COMPILE_PENCIL_FLAGS} -I${PENCIL_INCLUDE_DIRS} -o ${PENCIL_FILE_NAME}.ppcg.c ${CMAKE_CURRENT_SOURCE_DIR}/${pencil_file}
                       DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/${pencil_file}
                      )
    list(APPEND ${COMPILE_PENCIL_DEST_SOURCES} ${PENCIL_FILE_NAME}.ppcg.c)
  endforeach()
  set(${COMPILE_PENCIL_DEST_SOURCES}      ${${COMPILE_PENCIL_DEST_SOURCES}}      PARENT_SCOPE)
  set(${COMPILE_PENCIL_DEST_INCLUDE_DIRS} ${${COMPILE_PENCIL_DEST_INCLUDE_DIRS}} PARENT_SCOPE)
endfunction()
