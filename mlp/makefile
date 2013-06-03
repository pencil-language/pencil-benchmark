

all: test_mlp opencl_mlp

bench_mlp.gch: bench_mlp.hpp $(BASEDIR)/core/memory.gch
	$(CXX) $(CXXFLAGS) -c bench_mlp.hpp -o bench_mlp.gch

mlp_impl.o: mlp_impl.c

test_mlp.o: test_mlp.cpp bench_mlp.gch $(BASEDIR)/core/memory.gch

opencl_mlp.o: opencl_mlp.cpp bench_mlp.gch $(BASEDIR)/core/memory.gch

test_mlp: test_mlp.o mlp_impl.o allocator.o $(BASEDIR)/opencl/errors.o

opencl_mlp: opencl_mlp.o mlp_impl.o allocator.o  $(BASEDIR)/opencl/errors.o

allocator.gch: allocator.hpp $(BASEDIR)/core/memory.gch
	$(CXX) $(CXXFLAGS) -c allocator.hpp -o allocator.gch

allocator.o: allocator.cpp allocator.gch $(BASEDIR)/core/memory.gch

clean:
	-rm -f *.o test_mlp opencl_mlp *.gch test_operators
