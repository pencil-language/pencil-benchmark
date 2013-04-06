

all: test_mlp

mlp_impl.o: mlp_impl.c

test_mlp.o: test_mlp.cpp

test_mlp: test_mlp.o mlp_impl.o

clean:
	-rm -f *.o test_mlp
