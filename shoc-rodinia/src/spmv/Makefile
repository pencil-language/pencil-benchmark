include ../../config
include ../../scripts/default.mk

PPCG=$(PPCG_PATH)/ppcg

# Default PPCG flags, can be overriden in make invocation
PPCGFLAGS=--no-private-memory --no-shared-memory

# Benchmark invocation flags, can be overriden in make invocation
BENCHMARK_FLAGS=csr

CFLAGS+=-I$(PENCIL_RUNTIME) -I$(PENCIL_RUNTIME)/opencl-opt -I$(PPCG_PATH)

OBJS=spmv_host.o spmv-ref.o

ORIG_BIN=${SHOC_DIR}/bin/Serial/OpenCL/Spmv


all:

build: spmv-sp spmv-dp

ppcg: spmv.c
	$(PPCG) --opencl-pencil-runtime $(AUTOTUNER_PPCG_FLAGS) spmv.c

# Generate PENCIL version
spmv.c: gen_pencil.c spmv.h
	cpp -P $< -o $@

# Generate reference version
spmv-ref.c: gen_ref.c spmv.h
	cpp -P $< -o $@

spmv-ref.o: spmv-ref.c
	$(CC) -c $(CFLAGS) -include ${PENCIL_INCLUDE}/pencil.h $^ -o $@

# We need to compile util.c without -std=c99 to get drand48 and friends.
util-sp.o: util.c
	$(CC) -DDATATYPE=float -c $^ -o $@

util-dp.o: util.c
	$(CC) -DDATATYPE=double -c $^ -o $@

# Binaries for single and double precision
PENCIL_RT_FLAGS=-L$(PENCIL_RUNTIME)/lib -locl_pencil
spmv-sp: main.c util-sp.o $(OBJS) $(PPCG_PATH)/ocl_utilities.c
	$(CC) -DDATATYPE_SP $(CFLAGS) $(PENCIL_RT_FLAGS) $(LDFLAGS) -I ../../include $^ -o $@

spmv-dp: main.c util-dp.o $(OBJS) $(PPCG_PATH)/ocl_utilities.c
	$(CC) -DDATATYPE_DP $(CFLAGS) $(PENCIL_RT_FLAGS) $(LDFLAGS) -I ../../include $^ -o $@

UPLOAD_FILES=spmv-sp spmv-dp spmv_kernel.cl
UPLOAD_FILES_ORIG=$(ORIG_BIN) orig/run.sh

# Run targets

run-ellpackr: $(UPLOAD_FILES) ellpackr/run.sh
	../../scripts/run-on-board.sh $^

run-csr: $(UPLOAD_FILES) csr/run.sh
	../../scripts/run-on-board.sh $^

$(ORIG_BIN): .FORCE
	cp ../../include/measure-time.h ${SHOC_DIR}/src/opencl/level1/spmv
	${MAKE} -C ${SHOC_DIR}/src/common
	${MAKE} -C ${SHOC_DIR}/src/opencl/common
	$(MAKE) -C ${SHOC_DIR}/src/opencl/level1/spmv
	$(MAKE) -C ${SHOC_DIR}/src/opencl/level1/spmv install

clean:
	rm -f *.o spmv-sp spmv-dp spmv_host.c spmv_kernel.cl spmv_kernel.h spmv.c spmv-ref.c

include ../../scripts/benchmark-run.mk

run:
	@echo "Use make run-csr or make run-ellpackr"

.FORCE:
.PHONY: build ppcg clean all
