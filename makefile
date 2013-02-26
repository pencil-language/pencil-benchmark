
CFLAGS=-O0 -g

all: mlp_impl.o

mlp_impl.o: mlp_impl.c

clean:
	-rm -f *.o
