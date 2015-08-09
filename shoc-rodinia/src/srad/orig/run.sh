#!/bin/bash

set -e
set -o pipefail
./srad 100 0.5 502 458 | grep TOTAL_TIME | sed -e 's/TOTAL_TIME: //'
