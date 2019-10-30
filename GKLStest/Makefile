CXX=g++
CXXFLAGS=-c -O -std=c++11 -fopenmp -Wall
CC=gcc
CCFLAGS=-c

GKLStest: example.o rnd_gen.o gkls.o
	$(CXX) -o GKLStest -O example.o rnd_gen.o gkls.o -fopenmp

Rv: Rexample.o rnd_gen.o gkls.o
	$(CXX) -o R_GKLStest -O Rexample.o rnd_gen.o gkls.o -fopenmp

Rexample.o: Rexample.cpp
	$(CXX) $(CXXFLAGS) Rexample.cpp

example.o: example.cpp
	$(CXX) $(CXXFLAGS) example.cpp

gkls.o: gkls.c
	$(CC) $(CCFLAGS) gkls.c

rnd_gen.o: rnd_gen.c
	$(CC) $(CCFLAGS) rnd_gen.c

clean:
	rm -rf *.o GKLStest R_GKLStest
