CC=gcc
CXX=g++
CFLAGS=-O3 -Wall -lm -lgsl -lgslcblas
OBJ=angle.o file.o vector.o

%.o: %.c
	$(CC) -c -o $@ $< $(CFLAGS)

angle: $(OBJ)
	g++ -o $@ $^ $(CFLAGS)

clean:  
	rm -f angle *.o
