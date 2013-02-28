
CXXFLAGS= -O0 -g -I/home/ujoimro/Inst/opencv/build/OpenCV-2.4.2/optimized/install/include
CFLAGS=$(CXXFLAGS) -Wimplicit -Werror=implicit-function-declaration

LDFLAGS=-lboost_serialization -L/home/ujoimro/Inst/opencv/build/optimized/install/lib -lopencv_core

all: test_mlp

mlp_impl.o: mlp_impl.c

test_mlp: test_mlp.cpp mlp_impl.o cast.h mlp.hpp

clean:
	-rm -f *.o test_mlp
