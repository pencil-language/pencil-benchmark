include ../../config
include ../../scripts/default.mk

#C compiler
HOSTCC = gcc
C_SRC = main.c srad-pencil.c
C_FLAGS = -std=c99 -I ../../include
C_EXE = c_srad

PENCIL_SRC = srad-pencil.c
OCL_SRC = srad-pencil_host.c main.c $(PPCG_PATH)/ocl_utilities.c

ORIG_BENCH_DIR = ${RODINIA_DIR}/opencl/srad

ORIG_EXE = ${ORIG_BENCH_DIR}/srad
OCL_EXE = ocl_srad

ORIG_DATA = ${RODINIA_DIR}/data/srad/image.pgm

PENCIL_INCLUDE = ${PENCIL_UTIL_HOME}/include/pencil.h

pencil.h: ${PENCIL_INCLUDE}
	cpp -D__PENCIL__ -P $^ -o $@

build: $(OCL_EXE)

ppcg: $(PENCIL_SRC)
	$(PPCG) --opencl-pencil-runtime $(AUTOTUNER_PPCG_FLAGS) $^

$(PENCIL_SRC): srad-pencil.unopt.c pencil.h
	$(PENCIL_TOOLS_HOME)/optimizer $< -o $@ --include pencil.h

$(C_EXE): $(C_SRC)
	$(HOSTCC) $(C_FLAGS) $(C_SRC) -o $(C_EXE) -lm -DHOST_DEBUG

$(ORIG_EXE): .FORCE
	CPATH=${OPENCL_SDK}/include $(MAKE) -C ${ORIG_BENCH_DIR} C_C=${CC} OPENCL_LIB=${OPENCL_SDK}/lib CFLAGS="-I`pwd`/../../include"

%_host.c %_kernel.cl: $(PENCIL_SRC)
	$(PPCG) --target=opencl $(AUTOTUNER_PPCG_FLAGS) $^

$(OCL_EXE): srad-pencil_host.c
	$(CC) $(C_FLAGS) -I$(PPCG_PATH) -I$(PENCIL_RUNTIME) -I$(PENCIL_RUNTIME)/opencl-opt  -I$(OPENCL_SDK)/include/ $(OCL_SRC) -L$(OPENCL_SDK)/lib/ -lOpenCL -L$(PENCIL_RUNTIME)/lib -locl_pencil -lm -o $@

UPLOAD_FILES=$(OCL_EXE) srad-pencil_kernel.cl run.sh ${ORIG_DATA}
UPLOAD_FILES_ORIG=$(ORIG_EXE) ${ORIG_BENCH_DIR}/main.h ${ORIG_BENCH_DIR}/kernel/kernel_gpu_opencl.cl orig/run.sh ${ORIG_DATA}

run_c: $(C_EXE)
	./$(C_EXE) ${ORIG_DATA} 100 0.5 502 458

clean:
	rm -f $(C_EXE) $(OCL_EXE) srad-pencil_kernel.cl srad-pencil_kernel.h srad-pencil_host.c image_out.pgm ocl_image_out.pgm

include ../../scripts/benchmark-run.mk

.FORCE:

.PHONY: build ppcg clean run_c
