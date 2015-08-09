KERNEL=spmv

/home/b/src/ppcgs/pencilcc_prefix/bin/ppcg -D__PENCIL__ --target=prl ${KERNEL}.pencil.c -I/home/b/src/pencil_codes/pencil-headers/include/ -I /home/b/src/pencil_codes/prl_lib/include/ 

gcc --std=gnu99 -D__PENCIL__ test_${KERNEL}.c -I /home/b/src/pencil_codes/prl_lib/include/ -I/opt/AMDAPP/include/ -I../../include/ -Wl,-rpath=RIGIN:/home/b/src/ppcgs/pencilcc_prefix/lib/  -L/home/b/src/ppcgs/pencilcc_prefix/lib/ -lprl_opencl -L/opt/AMDAPP/lib/x86_64/ -lOpenCL -I/home/b/src/pencil_codes/pencil-headers/include/ ${KERNEL}.pencil_host.c -lm

echo "Executing Optimized Code"
./a.out csr
