CXX=g++
CXXFLAGS=-c -std=c++11 -fopenmp -Wall
CC=gcc
CCFLAGS=-c

GKLStest: example.o rnd_gen.o gkls.o
	$(CXX) -o GKLStest example.o rnd_gen.o gkls.o -fopenmp

example.o: example.cpp
	$(CXX) $(CXXFLAGS) example.cpp

gkls.o: gkls.c
	$(CC) $(CCFLAGS) gkls.c

rnd_gen.o: rnd_gen.c
	$(CC) $(CCFLAGS) rnd_gen.c

clean:
	rm -rf *.o GKLStest
