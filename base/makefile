# UjoImro, 2013

all: test_color expand test_nesting test_local test_allocation

test_color: test_color.o $(BASEDIR)/opencl/errors.o
	$(CXX) $(LDFLAGS) $(TARGET_ARCH) $^ $(LOADLIBES) $(LDLIBS) -o $@

test_color.o: test_color.cpp $(BASEDIR)/opencl/opencl.hpp $(BASEDIR)/core/utility.hpp 

test_nesting: test_nesting.o $(BASEDIR)/opencl/errors.o
	$(CXX) $(LDFLAGS) $(TARGET_ARCH) $^ $(LOADLIBES) $(LDLIBS) -o $@

test_nesting.o: test_nesting.cpp $(BASEDIR)/opencl/opencl.hpp $(BASEDIR)/core/utility.hpp

test_local: test_local.o $(BASEDIR)/opencl/errors.o
	$(CXX) $(LDFLAGS) $(TARGET_ARCH) $^ $(LOADLIBES) $(LDLIBS) -o $@

test_local.o: test_local.cpp $(BASEDIR)/opencl/opencl.hpp $(BASEDIR)/core/utility.hpp

test_allocation: test_allocation.o $(BASEDIR)/opencl/errors.o
	$(CXX) $(LDFLAGS) $(TARGET_ARCH) $^ $(LOADLIBES) $(LDLIBS) -o $@

test_allocation.o: test_allocation.cpp $(BASEDIR)/core/memory.hpp $(BASEDIR)/core/utility.hpp

expand: expand.o
	$(CXX) $(LDFLAGS) $(TARGET_ARCH) $^ $(LOADLIBES) $(LDLIBS) -o $@

expand.o: expand.cpp

clean: 
	-rm -f expand test_color *.o test_nesting test_local test_allocation

# LuM end of file
