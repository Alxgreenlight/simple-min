CXX=g++
CXXFLAGS=-c -std=c++11 -fopenmp
CC=gcc
CCFLAGS=-c

GKLStesttool: example_gkls.o rnd_gen.o gkls.o
	$(CXX) -o GKLStesttool example_gkls.o rnd_gen.o gkls.o -fopenmp

example_gkls.o: example_gkls.cpp
	$(CXX) $(CXXFLAGS) example_gkls.cpp

gkls.o: gkls.c
	$(CC) $(CCFLAGS) gkls.c

rnd_gen.o: rnd_gen.c
	$(CC) $(CCFLAGS) rnd_gen.c

clean:
	rm -rf *.o GKLStesttool