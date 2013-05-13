
all: memory.gch utility.gch

memory.gch: memory.hpp
	$(CXX) $(CXXFLAGS) -c memory.hpp -o memory.gch	

utility.gch: utility.hpp
	$(CXX) $(CXXFLAGS) -c utility.hpp -o utility.gch	

clean:
	-rm -f *.gch
