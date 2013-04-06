# UjoImro, 2013

CXXFLAGS=-std=c++0x

all: errors.gch opencl.gch

errors.gch: errors.hpp
	$(CXX) $(CXXFLAGS) -c errors.hpp -o errors.gch

opencl.gch: opencl.hpp
	$(CXX) $(CXXFLAGS) -c opencl.hpp -o opencl.gch

clean:
	-rm -f *.gch

# LuM end of file
