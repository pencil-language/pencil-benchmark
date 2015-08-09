#!/bin/bash

set -e
set -o pipefail
./ocl_srad image.pgm 100 0.5 502 458 | grep TOTAL_TIME | sed -e 's/TOTAL_TIME: //'
