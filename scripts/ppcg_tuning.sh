#!/bin/bash

# Set automatically a variable that points to the benchmark root (main)
# directory.
BENCHMARK_ROOT_DIRECTORY="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/../"

. $BENCHMARK_ROOT_DIRECTORY/scripts/scripts_config.conf

. $BENCHMARK_ROOT_DIRECTORY/scripts/core.sh
