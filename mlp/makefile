
CFLAGS=-O0 -g -I/home/ujoimro/Inst/opencv/build/OpenCV-2.4.2/optimized/install/include -Wimplicit -Werror=implicit-function-declaration -Werror=int-to-pointer-cast 
CXXFLAGS=$(CFLAGS) -std=c++0x

#LDFLAGS=-lboost_serialization -L/home/ujoimro/Inst/opencv/build/OpenCV-2.4.2/optimized/install/lib -lopencv_core -lefence
LDFLAGS=-lboost_serialization -L/home/ujoimro/Inst/opencv/build/OpenCV-2.4.2/optimized/install/lib -lopencv_core -lstdc++ -lm

all: test_mlp

mlp_impl.o: mlp_impl.c

test_mlp.o: test_mlp.cpp

test_mlp: test_mlp.o mlp_impl.o

clean:
	-rm -f *.o test_mlp
