#!/bin/bash

binary_path=$1
image=$2

BENCHMARK_ROOT_DIRECTORY="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/../"
BINARY_DIRECTORY=`dirname $binary_path`
BINARY_FILE=`basename $binary_path`


if [ -z $binary_path ]; then
	echo "./script <binary_path> <image_path>"
fi

if [ -z $image ]; then
	image=$BENCHMARK_ROOT_DIRECTORY/images/M104_ngc4594_sombrero_galaxy_hi-res.jpg
	echo "Default image images/M104_ngc4594_sombrero_galaxy_hi-res.jpg is being used"
fi

cd $BINARY_DIRECTORY

COMPUTE_PROFILE=2 COMPUTE_PROFILE_CSV=1 COMPUTE_PROFILE_CONFIG=$BENCHMARK_ROOT_DIRECTORY/scripts/Nvidia_profiling_options.conf COMPUTE_PROFILE_LOG=profiling_log_${BINARY_FILE}_%d.csv ./$BINARY_FILE $image

echo "The profiling information are in build/profiling_log_${BINARY_FILE}_*.csv"
echo "Due to the limited support of OpenCL in Nvidia drivers, only the profiling"
echo "of the first context in the program is possible (this means profiling for"
echo "the code generated by the PENCIL compiler)"
