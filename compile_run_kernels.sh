#!/bin/bash

#set -x

PPCG_COMPILER=~/src/ppcgs/ppcg-gforge-pet-tree/ppcg 
#PPCG_COMPILER=~/src/ppcgs/ppcg_summary/ppcg
BENCH_ROOT=~/src/pencil_codes/CARP-Benchmarks/
OPENCL_PREFIX=/usr/local/cuda-5.5/
#OPENCL_PREFIX=/opt/AMDAPP/

# Use PPCG pet tree for all except for filter2D and gaussian where PPCG summary must be used !

#ppcg_pet_tree:
PPCG_EXTRA_OPTIONS="--target=opencl  -D__PENCIL__  --opencl-print-time-measurements "

#ppcg_summary:
#PPCG_EXTRA_OPTIONS="--target=opencl  -D__PENCIL__  --opencl-print-kernels-time-measurements   --opencl-include-file=../pencil/pencil_math.h"

LIST_OF_KERNELS="resize dilate cvt_color warpAffine filter2D gaussian"


# Carotte
#PPCG_OPTIONS[0]="--isl-schedule-fuse=min                    --no-private-memory --sizes={kernel[i]->tile[128,128];kernel[i]->grid[464,53];kernel[i]->block[8,32]}"
#PPCG_OPTIONS[1]="--isl-schedule-fuse=min                    --no-private-memory --sizes={kernel[i]->tile[16,16];kernel[i]->grid[16,16];kernel[i]->block[16,16]}"
#PPCG_OPTIONS[2]="--isl-schedule-fuse=min --no-shared-memory --no-private-memory --sizes={kernel[i]->tile[32,32];kernel[i]->grid[928,27];kernel[i]->block[4,64]}"
#PPCG_OPTIONS[3]="--isl-schedule-fuse=min --no-shared-memory --no-private-memory --sizes={kernel[i]->tile[32,32];kernel[i]->grid[928,27];kernel[i]->block[4,64]}"
#PPCG_OPTIONS[4]="--isl-schedule-fuse=min                    --no-private-memory --sizes={kernel[i]->tile[32,8];kernel[i]->grid[116,53];kernel[i]->block[32,32]}"
#PPCG_OPTIONS[5]="--isl-schedule-fuse=min --no-shared-memory --no-private-memory --sizes={kernel[i]->tile[32,32];kernel[i]->grid[464,53];kernel[i]->block[8,32]}"


# Fermi
PPCG_OPTIONS[0]="--isl-schedule-fuse=min  --no-private-memory --sizes={kernel[i]->tile[64,64];kernel[i]->grid[717,26];kernel[i]->block[4,64]}"
PPCG_OPTIONS[1]="--isl-schedule-fuse=min  --no-private-memory --no-shared-memory --sizes={kernel[i]->tile[32,32];kernel[i]->grid[180,101];kernel[i]->block[16,16]}"
PPCG_OPTIONS[2]="--isl-schedule-fuse=min  --no-private-memory --no-shared-memory --sizes={kernel[i]->tile[64,64];kernel[i]->grid[180,101];kernel[i]->block[16,16]}"
PPCG_OPTIONS[3]="--isl-schedule-fuse=min  --no-private-memory --sizes={kernel[i]->tile[64,64];kernel[i]->grid[180,101];kernel[i]->block[16,16]}"
PPCG_OPTIONS[4]="--isl-schedule-fuse=min  --no-private-memory --no-shared-memory --sizes={kernel[i]->tile[32,32];kernel[i]->grid[717,26];kernel[i]->block[4,64]}"
PPCG_OPTIONS[5]="--isl-schedule-fuse=min  --no-private-memory --sizes={kernel[i]->tile[32,32];kernel[i]->grid[717,26];kernel[i]->block[4,64]}"


#LIST_OF_KERNELS="resize"
#PPCG_OPTIONS[0]="--isl-schedule-fuse=min  --no-private-memory  --sizes={kernel[i]->tile[64,64];kernel[i]->grid[717,26];kernel[i]->block[4,64]}"

NB_TESTS=3
PEROFRM_ONLY_ONE_TEST=1
COMPILE_WITH_PPCG_AND_USE_AUTOTUNING_OPTIONS=1
USE_DEFAULT_MAKE_FILE=0
LOG_FILE=log
######################################################################"
DEFINES=""
if [ $PEROFRM_ONLY_ONE_TEST = 1 ]; then
	DEFINES="$DEFINES -DRUN_ONLY_ONE_EXPERIMENT"
fi
######################################################################"

# INPUT: a boolean.
# Print success if input=0 (represents successful execution), else print error.
# USAGE: "success $?"
success()
{
  if [ $1 = 0 ]; then
    echo -e "\e[32m    .Success\e[0m"
  else
    echo -e "\e[31m    .Error\e[0m"
    ERROR_SOME_WHERE=1
  fi
}


# compile the kernel ($1) with ppcg and then with g++
compile()
{
  KERNEL=$1
  ppcg_options=$2
  echo
  echo "[$KERNEL]"
  echo "[$KERNEL]" >> $LOG_FILE

  if [ $COMPILE_WITH_PPCG_AND_USE_AUTOTUNING_OPTIONS = 1 ]; then
    echo "    .ppcg $ppcg_options"
    echo "    .ppcg $PPCG_EXTRA_OPTIONS $ppcg_options $KERNEL.pencil.c" >> $LOG_FILE
    $PPCG_COMPILER $PPCG_EXTRA_OPTIONS $ppcg_options -I$BENCH_ROOT/$KERNEL $BENCH_ROOT/$KERNEL/$KERNEL.pencil.c &>> $LOG_FILE
    success $?
  fi

  echo "    .compiling ${KERNEL}.pencil_host.c and test_${KERNEL}.cpp (g++)"
  g++ -x c -c -O3 -DNDEBUG -march=native -fomit-frame-pointer -fPIC -ffast-math -Wall -std=c99 -Iinclude -Ibuild -I$OPENCL_PREFIX/include/ -I$BENCH_ROOT/$KERNEL ${KERNEL}.pencil_host.c -o $KERNEL.pencil_host.o &>> $LOG_FILE
  EXIT_STATUS_COMPILATION_1=$?

  g++ -shared -O3 -o lib${KERNEL}_ppcg.so $BENCH_ROOT/build/ocl_utilities.o $KERNEL.pencil_host.o -L$OPENCL_PREFIX/lib/x86_64/ -lOpenCL -L/usr/local/lib/ -lopencv_core -lopencv_imgproc -lopencv_ocl -lopencv_highgui -L/usr/lib/ -ltbb -ltbbmalloc -Lbuild -lboost_date_time -lboost_filesystem -lboost_iostreams -lboost_program_options -lboost_serialization -lboost_system -lboost_chrono &>> $LOG_FILE
  EXIT_STATUS_COMPILATION_2=$?

  g++ -O3 $DEFINES -DNDEBUG -march=native -fomit-frame-pointer -fPIC -ffast-math -Wall -std=c++0x -Iinclude -Ibuild -I$OPENCL_PREFIX/include/ -I/usr/local/include/ -I/usr/include/ -I$BENCH_ROOT/include/ -Wl,-rpath=RIGIN:/usr/local/lib/ -I$BENCH_ROOT/build/ $BENCH_ROOT/$KERNEL/test_${KERNEL}.cpp -o ppcg_test_${KERNEL} -L$OPENCL_PREFIX/lib/x86_64/ -lOpenCL -L/usr/local/lib/ -lopencv_core -lopencv_imgproc -lopencv_ocl -lopencv_highgui -L/usr/lib/ -ltbb -ltbbmalloc -Lbuild -lboost_date_time -lboost_filesystem -lboost_iostreams -lboost_program_options -lboost_serialization -lboost_system -lboost_chrono -L. -l${KERNEL}_ppcg &>> $LOG_FILE
  EXIT_STATUS_COMPILATION_3=$?

  EXIT_STATUS_COMPILATION=`expr $EXIT_STATUS_COMPILATION_1 + $EXIT_STATUS_COMPILATION_2 + $EXIT_STATUS_COMPILATION_3`
  success $EXIT_STATUS_COMPILATION
}


run()
{
  KERNEL=$1

  echo -n "    .running ./ppcg_test_${KERNEL}: "
 
  for ((i=0; i < $NB_TESTS; i++)); do
	  echo -n "$i/$NB_TESTS "
	  ./ppcg_test_${KERNEL} 
	  ppcg_exit_status=$?
  done
  echo

  success $ppcg_exit_status

  echo "--------------------------------------------------" >> $LOG_FILE
}

####################################################################################"

if [ $USE_DEFAULT_MAKE_FILE = 1 ]; then
	make clean
	make -j4
	clear
	echo "Running benchmarks"
	echo "-------------------------------------------------------"
fi

./prepare_pool.sh
cd build
rm -rf $LOG_FILE

id=0;
for ker in ${LIST_OF_KERNELS}; do
	options="$PPCG_EXTRA_OPTIONS ${PPCG_OPTIONS[$id]}"
        compile $ker "$options";
	id=`expr $id + 1`
        run $ker;
done

cd ..
./restore_pool.sh
