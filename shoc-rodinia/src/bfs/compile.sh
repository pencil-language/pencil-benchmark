KERNEL=bfs

gcc --std=c99 test_${KERNEL}.c -I/home/b/src/ppcgs/pencilcc_prefix/include/ -I/opt/AMDAPP/include/ -I../../include/ -Wl,-rpath=RIGIN:/home/b/src/ppcgs/pencilcc_prefix/lib/ -L/home/b/src/ppcgs/pencilcc_prefix/lib/ -lprl_opencl -L/opt/AMDAPP/lib/x86_64/ -lOpenCL -I/home/b/src/pencil_codes/pencil-headers/include/  ${KERNEL}.pencil.c  -lm

echo "Executing native"
./a.out graph4M.txt

/home/b/src/ppcgs/pencilcc_prefix/bin/ppcg -D__PENCIL__ --target=prl ${KERNEL}.pencil.c -I/home/b/src/pencil_codes/pencil-headers/include/ -I/home/b/src/ppcgs/pencilcc_prefix/include/
gcc --std=c99 test_${KERNEL}.c -I/home/b/src/ppcgs/pencilcc_prefix/include/ -I/opt/AMDAPP/include/ -I../../include/ -Wl,-rpath=RIGIN:/home/b/src/ppcgs/pencilcc_prefix/lib/  -L/home/b/src/ppcgs/pencilcc_prefix/lib/ -lprl_opencl -L/opt/AMDAPP/lib/x86_64/ -lOpenCL ${KERNEL}.pencil_host.c -I/home/b/src/pencil_codes/pencil-headers/include/ -lm

echo "Executing Optimized Code"
./a.out graph4M.txt

