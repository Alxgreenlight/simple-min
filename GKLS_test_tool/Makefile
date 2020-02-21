CXX=g++
CXXFLAGS=-c -std=c++11 -O -Wall -fopenmp
CC=gcc
CCFLAGS=-c
INCDIR=-I../
GKLSDIR=../GKLStest/

GKLStesttool: example_gkls.o rnd_gen.o gkls.o
	$(CXX) -o GKLStesttool example_gkls.o rnd_gen.o gkls.o -fopenmp

example_gkls.o: example_gkls.cpp
	$(CXX) $(INCDIR) $(CXXFLAGS) example_gkls.cpp

gkls.o:
	$(CC) $(CCFLAGS) "$(GKLSDIR)gkls.c"

rnd_gen.o:
	$(CC) $(CCFLAGS) "$(GKLSDIR)rnd_gen.c"

clean:
	rm -rf *.o GKLStesttool

rebuild: clean GKLStesttool