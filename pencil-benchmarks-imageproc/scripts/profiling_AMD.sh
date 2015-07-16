#!/bin/bash

# Example PROFILER=/opt/AMD/CodeXL_1.6-7247/x86_64/sprofile
PROFILER=

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

$PROFILER -p -O -o profiling_AMD_log_${BINARY_FILE}.csv ${BINARY_FILE}  $image

echo "The profiling information are in build/profiling_AMD_log_${BINARY_FILE}.csv"


