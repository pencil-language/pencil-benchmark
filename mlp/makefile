

all: test_mlp opencl_mlp

bench_mlp.gch: bench_mlp.hpp
	$(CXX) $(CXXFLAGS) -c bench_mlp.hpp -o bench_mlp.gch

mlp_impl.o: mlp_impl.c

test_mlp.o: test_mlp.cpp bench_mlp.gch

opencl_mlp.o: opencl_mlp.cpp bench_mlp.gch

test_mlp: test_mlp.o mlp_impl.o 

opencl_mlp: opencl_mlp.o mlp_impl.o

clean:
	-rm -f *.o test_mlp opencl_mlp *.gch
