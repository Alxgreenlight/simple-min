CXX=g++
CXXFLAGS=-c -std=c++11 -fopenmp -Wall
CC=gcc
CCFLAGS=-c

miniOneDim: OneDimStarter.o
	$(CXX) -o miniOneDim OneDimStarter.o -fopenmp

OneDimStarter.o: OneDimStarter.cpp
	$(CXX) $(CXXFLAGS) OneDimStarter.cpp

clean:
	rm -rf *.o miniOneDim