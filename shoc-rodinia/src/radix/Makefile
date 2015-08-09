BENCH_SRC = ${SHOC_DIR}/src/opencl/level1/sort/
BIN = radix
OCL_UTILITIES = ${PPCG_PATH}/ocl_utilities.c
N_NUMBERS = 16777216
ORIG_SIZE = 4
TEST_FILE = n${N_NUMBERS}
ORIG_BIN = ${SHOC_DIR}/bin/Serial/OpenCL/Sort

UPLOAD_FILES = pencil_kernel.cl ${BIN} run.sh
UPLOAD_FILES_ORIG = ${ORIG_BIN} orig/run.sh

DATA_FILES = ${TEST_FILE}

include ../../config
include ../../scripts/default.mk
include ../../scripts/benchmark-run.mk

ppcg:pencil.c
	${PPCG_PATH}/ppcg ${AUTOTUNER_PPCG_FLAGS} --opencl-compiler-options="-Dassume(x)" $^

radix: pencil_host.c ${OCL_UTILITIES} sort.c
	${CC} ${CFLAGS} -I$(PPCG_PATH) ${LDFLAGS} $^ -o $@ -I. -I../../include -include ${PENCIL_INCLUDE}/pencil.h -I$(PENCIL_RUNTIME) -I$(PENCIL_RUNTIME)/opencl-opt -L$(PENCIL_RUNTIME)/lib -locl_pencil

gen_tests: gen_tests.c
	gcc -std=c99 $^ -o $@

${TEST_FILE}: gen_tests
	./gen_tests ${N_NUMBERS}

${ORIG_BIN}: .FORCE
	make -C ${SHOC_DIR}/src/common
	make -C ${SHOC_DIR}/src/opencl/common
	make -C ${BENCH_SRC} clean
	make -C ${BENCH_SRC}
	make -C ${BENCH_SRC} install

build: radix ${TEST_FILE}

pencil.unopt.c: pencil.base.c
	cpp -P $^ > $@

pencil.c: pencil.unopt.c
	${PENCIL_TOOLS_HOME}/optimizer $^ -o $@

init: pencil.c

.FORCE:

.PHONY: run.sh
run.sh:
	echo "#!/bin/bash" > $@
	echo "set -e" >> $@
	echo "set -o pipefail" >> $@
	echo "./radix ${BOARD_DATA_DIR}/${TEST_FILE} | grep TOTAL_TIME | sed -e 's/TOTAL_TIME: //'" >> $@
	chmod a+x $@

.PHONY: orig/run.sh
orig/run.sh:
	echo "#!/bin/bash" > $@
	echo "set -e" >> $@
	echo "set -o pipefail" >> $@
	echo "./Sort --passes 1 --size ${ORIG_SIZE} | grep TOTAL_TIME | sed -e 's/TOTAL_TIME: //'" >> $@
	chmod a+x $@
