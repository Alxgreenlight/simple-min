CXX=g++
CXXFLAGS=-c -std=c++11 -fopenmp -Wall
CC=gcc
CCFLAGS=-c

mathexptest: tutorial.o
	$(CXX) -o mathexptest tutorial.o -fopenmp

tutorial.o: tutorial.cpp
	$(CXX) $(CXXFLAGS) tutorial.cpp -fopenmp

clean:
	rm -rf *.o mathexptest