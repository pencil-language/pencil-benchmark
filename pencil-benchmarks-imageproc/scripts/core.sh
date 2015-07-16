#!/bin/bash

# This script should not be called directly.  It should be called using
# other wrappers.
# Use this script to compile the benchmark kernels and to perform tuning.
# To perform tuning run the script without any option.
# In this case, the script will generate many PPCG options, and will
# launche PPCG with each one of these options, then will compile each
# output of PPCG and will execute it.
# To compile the kernels using preset PPCG options (obtained from
# previous tuning runs):
# ./compile_and_tune_ppcg.sh <ppcg_options_file>
# The folder scripts/ppcg_preset_options/ contains files with preset
# PPCG options that were obtained through tuning on different
# architectures.
# <ppcg_options_file> contains the list of optimization options
# for the kernels.  The option that corresponds to the i'th kernel
# in $LIST_OF_KERNELS (a variable defined in this script) is stored
# in best_optimization_options[$i].
# The execution times for each PPCG options for
# a given kernel are reported in build/output_time.$KERNEL_NAME.csv.
# The best execution times for each kernel are reported in
# build/$OUTPUT_TIME_FILE.csv
# The best PPCG options obtained through tuning are reported in
# the file build/best_optimizations_log.sh
# The log for kernels compilation and execution is stored in build/$LOG_FILE


# Perform only one experiment for each kernel (i.e. run one iteration
# and use one configuration only).
DEFINED_VARIABLES="-DRUN_ONLY_ONE_EXPERIMENT"

if [ -z "$1" ]; then
	BEST_OPTIMIZATIONS_LOG=$BENCHMARK_ROOT_DIRECTORY/build/best_optimizations_log.sh
	ENABLE_TUNING=1
else
	BEST_OPTIMIZATIONS_LOG=$1
	ENABLE_TUNING=0
fi

OUTPUT_TIME_FILE="output_time"
TEMP_OUTPUT_FILE=temp_output_file
TEMP_TIME_FILE_1=temp_time_file_1
TEMP_TIME_FILE_2=temp_time_file_2
TEMP_TIME_FILE_3=temp_time_file_3
LOG_FILE=benchmark_building_log.txt
CSV_DELIMITER="/"

if [ $USE_ARM_MALI_OPENCL_LIBRARIES = 0 ]; then
	OPENCL_LIBRARY="-lOpenCL"
else
	OPENCL_LIBRARY="-lmali"
fi

HEADER_FLAGS="-I$PENCIL_INCLUDE_DIR/ -I$PRL_INCLUDE_DIR/ -I$OPENCL_INCLUDE_DIR/ -I$OPENCV_INCLUDE_DIR/ -I$BENCHMARK_ROOT_DIRECTORY/include/"
LINKER_FLAGS="-L$PRL_LIB_DIR/ -L$OPENCL_LIB_DIR/ -L$BENCHMARK_ROOT_DIRECTORY/build/ -L$OPENCV_LIB_DIR/"
LIBRARY_FLAGS="-lprl -lopencv_core -lopencv_imgproc -lopencv_ocl -lopencv_highgui $OPENCL_LIBRARY -ltbb -ltbbmalloc"

PENCIL_COMPILER_EXTRA_OPTIONS="--target=opencl -D__PENCIL__ --opencl-include-file=${PENCIL_INCLUDE_DIR}/pencil_opencl.h"

if [ $TUNE_LOOP_FUSION_HEURISTICS = 1 ]; then
	POSSIBLE_LOOP_FUSION_OPTIONS[0]="--isl-schedule-fuse=max --no-isl-schedule-separate-components"
	POSSIBLE_LOOP_FUSION_OPTIONS[1]="--isl-schedule-fuse=min"
else
	POSSIBLE_LOOP_FUSION_OPTIONS[0]="--isl-schedule-fuse=min"
fi

if [ $TUNE_SHARED_MEMORY = 1 ]; then
	POSSIBLE_SHARED_MEM_OPTIONS[0]=""
	POSSIBLE_SHARED_MEM_OPTIONS[1]="--no-shared-memory"
else
	POSSIBLE_SHARED_MEM_OPTIONS[0]="--no-shared-memory"
fi

if [ $TUNE_PRIVATE_MEMORY = 1 ]; then
	POSSIBLE_PRIVATE_MEM_OPTIONS[0]=""
	POSSIBLE_PRIVATE_MEM_OPTIONS[1]="--no-private-memory"
else
	POSSIBLE_PRIVATE_MEM_OPTIONS[0]="--no-private-memory"
fi

###############################################################################

# Print Success (green) if the input ($1) is 0 (0 represents a successful execution),
# otherwise print error (red).
# USAGE: "test_success $?"
test_success()
{
  if [ $1 = 0 ]; then
    echo -e "\e[32m    .Success\e[0m"
  else
    echo -e "\e[31m    .Error\e[0m"
  fi
}

# Get the median of the list of numbers listed in the input file ($1).
# The numbers in the input file should be separated with a new line.
get_median()
{
	FILE=$1
	sort -n $FILE | nawk 'NF{a[NR]=$1;c++}END {printf (c%2==0)?(a[int(c/2)+1]+a[int(c/2)])/2:a[int(c/2)+1]}'
}

# Compile the kernel ($1) with ppcg and then with g++
compile()
{
  KERNEL=$1
  ppcg_tuning_options=$2
  echo
  echo "[$KERNEL]"
  echo "[$KERNEL]" >> $LOG_FILE

  echo "    .ppcg $ppcg_tuning_options"
  echo "    .ppcg $PENCIL_COMPILER_EXTRA_OPTIONS $ppcg_tuning_options $KERNEL.pencil.c" >> $LOG_FILE
  $PENCIL_COMPILER_BINARY $PENCIL_COMPILER_EXTRA_OPTIONS $ppcg_tuning_options $HEADER_FLAGS -I$BENCHMARK_ROOT_DIRECTORY/$KERNEL $BENCHMARK_ROOT_DIRECTORY/$KERNEL/$KERNEL.pencil.c &>> $LOG_FILE
  test_success $?

  echo "    .compiling ${KERNEL}.pencil_host.c and test_${KERNEL}.cpp (g++)"
  g++ -x c -c -O3 -DNDEBUG -fomit-frame-pointer -fPIC -std=c99 $HEADER_FLAGS -I$BENCHMARK_ROOT_DIRECTORY/$KERNEL ${KERNEL}.pencil_host.c -o $KERNEL.pencil_host.o &>> $LOG_FILE
  EXIT_STATUS_COMPILATION_1=$?

  g++ -shared -O3 -o lib${KERNEL}_ppcg.so $KERNEL.pencil_host.o $LINKER_FLAGS $LIBRARY_FLAGS &>> $LOG_FILE
  EXIT_STATUS_COMPILATION_2=$?

  g++ -O3 $DEFINED_VARIABLES -fomit-frame-pointer -fPIC -std=c++0x $HEADER_FLAGS -Wl,-rpath=RIGIN:$PRL_LIB_DIR $BENCHMARK_ROOT_DIRECTORY/$KERNEL/test_${KERNEL}.cpp -o ppcg_test_${KERNEL} $LIBRARY_FLAGS $LINKER_FLAGS -l${KERNEL}_ppcg &>> $LOG_FILE
  EXIT_STATUS_COMPILATION_3=$?

  EXIT_STATUS_COMPILATION=`expr $EXIT_STATUS_COMPILATION_1 + $EXIT_STATUS_COMPILATION_2 + $EXIT_STATUS_COMPILATION_3`
  test_success $EXIT_STATUS_COMPILATION
}

# Global variables used to keep track of the best optimization time
# and best optimization options.
BEST_EXECUTION_TIME=9999
BEST_OPTIMIZATION_OPTIONS=""

run()
{
  KERNEL=$1
  OPTIONS=$2
  OUTPUT_FILE=$3

  rm -rf $TEMP_TIME_FILE_1 $TEMP_TIME_FILE_2 $TEMP_TIME_FILE_3 $TEMP_OUTPUT_FILE

  echo -n "    .running ./ppcg_test_${KERNEL}: "

  for ((i=0; i < $NB_RUNS; i++)); do
	  echo -n "$i/$NB_RUNS "
	  PRL_BLOCKING=1 ./ppcg_test_${KERNEL} $TEST_IMAGE 1>>$TEMP_OUTPUT_FILE 2>>$LOG_FILE
	  exit_status=$?
  done
  echo

  test_success $exit_status
  if [ $exit_status = 0 ]; then

      # OpenCV total execution time without compilation time
      cat $TEMP_OUTPUT_FILE | grep -F "[RealEyes] Accumulate GPU time (inc copy):" | awk '{print $NF;}' 1> $TEMP_TIME_FILE_1
	  echo -n `get_median $TEMP_TIME_FILE_1` >> $OUTPUT_FILE
	  echo -n "$CSV_DELIMITER" >> $OUTPUT_FILE

      # Total execution time (computed from host) without compilation time
      cat $TEMP_OUTPUT_FILE | grep -F "Duration:" | awk '{print $2;}' 1> $TEMP_TIME_FILE_2
	  echo -n `get_median $TEMP_TIME_FILE_2` >> $OUTPUT_FILE
	  echo -n "$CSV_DELIMITER" >> $OUTPUT_FILE
	  current_execution_time=`get_median $TEMP_TIME_FILE_2`
	  if [ $(echo " $current_execution_time < $BEST_EXECUTION_TIME" | bc) -eq 1 ]; then
		BEST_EXECUTION_TIME=$current_execution_time
		BEST_OPTIMIZATION_OPTIONS=$OPTIONS
	  fi

      # Kernel only execution time computed using OpenCL profiling
      cat $TEMP_OUTPUT_FILE | grep -F "GPU_Compute:" | awk '{print $2;}' 1> $TEMP_TIME_FILE_3
	  echo `get_median $TEMP_TIME_FILE_3` >>  $OUTPUT_FILE
  else
	  echo " ERROR in ./ppcg_test_${KERNEL}" >> $LOG_FILE
	  echo "9999 $CSV_DELIMITER 9999 $CSV_DELIMITER 9999" >> $OUTPUT_FILE
  fi

  echo "--------------------------------------------------" >> $LOG_FILE
}

####################################################################################"

PREPARE_KERNEL_SPECIFIC_OUTPUT_FILE()
{
	KERNEL=$1

	rm -rf ${OUTPUT_TIME_FILE}.${KERNEL}.csv

	echo -n "PPCG options $CSV_DELIMITER" > ${OUTPUT_TIME_FILE}.${KERNEL}.csv

	echo -n "total_execution_time_opencv $CSV_DELIMITER" >> ${OUTPUT_TIME_FILE}.${KERNEL}.csv
	echo -n "total_execution_time_ppcg $CSV_DELIMITER" >> ${OUTPUT_TIME_FILE}.${KERNEL}.csv
	echo -n "kernel_only_execution_time_ppcg $CSV_DELIMITER" >> ${OUTPUT_TIME_FILE}.${KERNEL}.csv

	echo "" >> ${OUTPUT_TIME_FILE}.${KERNEL}.csv
}


PREPARE_GENERAL_OUTPUT_FILE()
{
	rm -rf ${OUTPUT_TIME_FILE}.csv

	echo "kernel $CSV_DELIMITER total_execution_time_opencv $CSV_DELIMITER total_execution_time_ppcg $CSV_DELIMITER kernel_only_execution_time_ppcg" >> ${OUTPUT_TIME_FILE}.csv

}

##########################

KERNEL_ID=0

dump_best_optimization_options()
{
   KERNEL=$1
   best_optimization_options=$2


   echo "best_optimization_options[$KERNEL_ID]=\"$best_optimization_options\"" >> $BEST_OPTIMIZATIONS_LOG
   KERNEL_ID=`expr $KERNEL_ID + 1`
}


AUTO_TUNE()
{
	KERNEL=$1

	BEST_EXECUTION_TIME=9999
	BEST_OPTIMIZATION_OPTIONS=""

	NB_TEST_1=`echo ${POSSIBLE_BLOCK_SIZES} | wc -w`
	NB_TEST_2=`echo ${POSSIBLE_TILE_SIZES}| wc -w`

	if [ $TUNE_TILE_AND_BLOCK_SIZES = 1 ]; then
		TOTAL_NUMBER_OF_TESTS=`expr ${#POSSIBLE_LOOP_FUSION_OPTIONS[@]} \* ${#POSSIBLE_SHARED_MEM_OPTIONS[@]} \* ${#POSSIBLE_PRIVATE_MEM_OPTIONS[@]} + ${#POSSIBLE_LOOP_FUSION_OPTIONS[@]} \* ${#POSSIBLE_SHARED_MEM_OPTIONS[@]} \* ${#POSSIBLE_PRIVATE_MEM_OPTIONS[@]} \* $NB_TEST_1 \* $NB_TEST_2`
	else
		TOTAL_NUMBER_OF_TESTS=`expr ${#POSSIBLE_LOOP_FUSION_OPTIONS[@]} \* ${#POSSIBLE_SHARED_MEM_OPTIONS[@]} \* ${#POSSIBLE_PRIVATE_MEM_OPTIONS[@]}`
	fi

	option_0=""
	option_1=""
	option_2=""
	option_3=""
	option_4=""
	option_counter=0

	  for i0 in ${!POSSIBLE_LOOP_FUSION_OPTIONS[*]}; do
	    for i1 in ${!POSSIBLE_SHARED_MEM_OPTIONS[*]}; do
	      for i2 in ${!POSSIBLE_PRIVATE_MEM_OPTIONS[*]}; do

		      option_counter=`expr $option_counter + 1`
		      option_0=${POSSIBLE_LOOP_FUSION_OPTIONS[$i0]}
		      option_1=${POSSIBLE_SHARED_MEM_OPTIONS[$i1]}
		      option_2=${POSSIBLE_PRIVATE_MEM_OPTIONS[$i2]}
		      options="$option_0 $option_1 $option_2"
		      echo
		      echo "Options [$option_counter/$TOTAL_NUMBER_OF_TESTS]: $options"

		      echo -n "$options $CSV_DELIMITER" >> ${OUTPUT_TIME_FILE}.${KERNEL}.csv
		      compile $KERNEL "$options";
		      run $KERNEL "$options" ${OUTPUT_TIME_FILE}.${KERNEL}.csv;

		      echo "" >> ${OUTPUT_TIME_FILE}.${KERNEL}.csv

		      if [ $TUNE_TILE_AND_BLOCK_SIZES = 1 ]; then
		        block_indice=1
			# Generate the different combinations for --sizes '{...}'
			for i5 in ${POSSIBLE_BLOCK_SIZES}; do
			  for i3 in ${POSSIBLE_TILE_SIZES}; do
		              option_counter=`expr $option_counter + 1`
			      option_0=${POSSIBLE_LOOP_FUSION_OPTIONS[$i0]}
			      option_1=${POSSIBLE_SHARED_MEM_OPTIONS[$i1]}
			      option_2=${POSSIBLE_PRIVATE_MEM_OPTIONS[$i2]}
			      option_4=${TUNING_OPENCL_COMPILER_OPTIONS[$i6]}
			      option_3="--sizes={kernel[i]->tile[$i3];kernel[i]->block[$i5]}"
			      options="$option_0 $option_1 $option_2 $option_3 $option_4"
			      echo
			      echo "Options [$option_counter/$TOTAL_NUMBER_OF_TESTS]: $options"
			      echo -n "$options $CSV_DELIMITER" >> ${OUTPUT_TIME_FILE}.${KERNEL}.csv

			      compile $KERNEL "$options";
		              run $KERNEL "$options" ${OUTPUT_TIME_FILE}.${KERNEL}.csv;

			      echo "" >> ${OUTPUT_TIME_FILE}.${KERNEL}.csv
			  done
			  block_indice=`expr $block_indice + 1`;
			done
		      fi # $TUNE_TILE_AND_BLOCK_SIZES = 1
	      done
	    done
	  done

	dump_best_optimization_options $KERNEL "$BEST_OPTIMIZATION_OPTIONS"
	echo "--------------------------------------------------"
}

# Setup LD_LIBRARY_PATH
LD_LIBRARY_PATH=.:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH

echo
echo "Timings will be generated in build/${OUTPUT_TIME_FILE}.csv"

if [ ! -d $BENCHMARK_ROOT_DIRECTORY/build ]; then
       mkdir $BENCHMARK_ROOT_DIRECTORY/build
fi

# Load best PPCG optimizations passed by the user.
if  [ $ENABLE_TUNING = 0 ]; then
	if [ -f $BEST_OPTIMIZATIONS_LOG ]; then
		source $BEST_OPTIMIZATIONS_LOG
	fi
fi

cd $BENCHMARK_ROOT_DIRECTORY/build
rm -rf $LOG_FILE

# Do autotuning if ENABLE_TUNING == 1.
if [ $ENABLE_TUNING = 1 ]; then
	rm -rf $BEST_OPTIMIZATIONS_LOG
	for ker in ${LIST_OF_KERNELS}; do
		PREPARE_KERNEL_SPECIFIC_OUTPUT_FILE $ker;
		AUTO_TUNE $ker
	done
	echo "Autotuning DONE" >> $LOG_FILE
	echo "****************" >> $LOG_FILE
	# Load the best PPCG options found by autotuning
	source $BEST_OPTIMIZATIONS_LOG
fi

PREPARE_GENERAL_OUTPUT_FILE;

# Copy the hog.pencil.cl file into the build directory
if [ ! -d "$BENCHMARK_ROOT_DIRECTORY/build/hog" ]; then
	mkdir $BENCHMARK_ROOT_DIRECTORY/build/hog
fi
cp $BENCHMARK_ROOT_DIRECTORY/hog/hog.opencl.cl $BENCHMARK_ROOT_DIRECTORY/build/hog

id=0
for ker in ${LIST_OF_KERNELS}; do
	echo -n "$ker $CSV_DELIMITER" >> ${OUTPUT_TIME_FILE}.csv
	options="${best_optimization_options[$id]}"
	compile $ker "$options"
	run $ker "$options" ${OUTPUT_TIME_FILE}.csv
	id=`expr $id + 1`
done

echo "DONE" >> $LOG_FILE
echo "****************" >> $LOG_FILE

cd $BENCHMARK_ROOT_DIRECTORY
