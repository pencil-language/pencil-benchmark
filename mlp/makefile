

TO_BUILD=test_mlp opencl_mlp

ifneq ($(NO_INTEL_LIBS), 1)
	TO_BUILD+=gel_mlp GEL
endif

all: $(TO_BUILD)

GEL:
	+make GEL

CFLAGS+=-std=gnu99

#bench_mlp.gch: bench_mlp.hpp $(BASEDIR)/core/memory.gch
#	$(CXX) $(CXXFLAGS) -c bench_mlp.hpp -o bench_mlp.gch

mlp_impl.o: mlp_impl.c

test_mlp.o: test_mlp.cpp bench_mlp.hpp $(BASEDIR)/core/memory.hpp

opencl_mlp.o: opencl_mlp.cpp bench_mlp.hpp $(BASEDIR)/core/memory.hpp

gel_mlp.o: gel_mlp.cpp GEL/histo.h GEL/histo_impl.h GEL/linalg.h GEL/defines.h

test_mlp: test_mlp.o mlp_impl.o allocator.o $(BASEDIR)/opencl/errors.o

opencl_mlp: opencl_mlp.o mlp_impl.o allocator.o  $(BASEDIR)/opencl/errors.o

gel_mlp:  gel_mlp.o mlp_impl.o allocator.o $(BASEDIR)/opencl/errors.o GEL/linalg.o 

#allocator.gch: allocator.hpp $(BASEDIR)/core/memory.gch
#	$(CXX) $(CXXFLAGS) -c allocator.hpp -o allocator.gch

allocator.o: allocator.cpp allocator.hpp $(BASEDIR)/core/memory.hpp

clean:
	-rm -f *.o test_mlp opencl_mlp gel_mlp *.gch test_operators
	make -C GEL clean



