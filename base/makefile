# UjoImro, 2013


all: test_color expand 


test_color: test_color.o 

test_color.o: test_color.cpp $(BASEDIR)/opencl/opencl.gch

expand: expand.o

expand.o: expand.cpp

clean: 
	-rm -f expand color *.o test_color

# LuM end of file
