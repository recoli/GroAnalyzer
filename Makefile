CC=gcc
CXX=g++
CFLAGS=-O3 -Wall -lm -lgsl -lgslcblas

%.o: %.c
	$(CC) -c -o $@ $< $(CFLAGS)

angle: angle.o file.o vector.o
	g++ -o $@ $^ $(CFLAGS)

op: op.o file.o vector.o
	g++ -o $@ $^ $(CFLAGS)

rotate: rotate.o elg.o file.o vector.o
	g++ -o $@ $^ $(CFLAGS)

clean:  
	rm -f angle op rotate *.o
