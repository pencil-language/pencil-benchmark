
CXXFLAGS= -O0 -g -I/home/ujoimro/Inst/opencv/build/OpenCV-2.4.2/optimized/install/include
CFLAGS=$(CXXFLAGS) -Wimplicit -Werror=implicit-function-declaration


all: test_mlp

mlp_impl.o: mlp_impl.c

test_mlp: test_mlp.cpp mlp_impl.o

clean:
	-rm -f *.o test_mlp
