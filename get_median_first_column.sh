#!/bin/bash

if [ $# -ne 1 ]; then
    echo "Usage: ./time_benchmarh.sh <file.txt>";
    echo "See example in the end of the script";
    exit 1;
fi;


sort -n $1 | nawk 'NF{a[NR]=$1;c++}END {printf (c%2==0)?(a[int(c/2)+1]+a[int(c/2)])/2:a[int(c/2)+1]}'


#example
#12.5
#13.6
#11.3
#16.2
