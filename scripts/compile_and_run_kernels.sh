#!/bin/bash

ppcg_options=$1

# Set automatically a variable that points to the benchmark root (main)
# directory.
BENCHMARK_ROOT_DIRECTORY="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/../"

. $BENCHMARK_ROOT_DIRECTORY/scripts/scripts_config.conf

echo

if [ -z $ppcg_options ]; then
	OPTIONS=$BENCHMARK_ROOT_DIRECTORY/scripts/ppcg_preset_options/ppcg_default_options.sh
	echo "Using PPCG default options."
else
	OPTIONS=$ppcg_options
	echo "Using PPCG preset options: $ppcg_options"
fi

. $BENCHMARK_ROOT_DIRECTORY/scripts/core.sh $OPTIONS
