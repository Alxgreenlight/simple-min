CXX=g++
CXXFLAGS=-c -std=c++11 -Wall -fopenmp
CC=gcc
CCFLAGS=-c

GrishaginTest: example.o vagris.o
	$(CXX) -o GrishaginTest example.o vagris.o -fopenmp

example.o: example.cpp
	$(CXX) $(CXXFLAGS) example.cpp

vagris.o: vagris.c
	$(CC) $(CCFLAGS) vagris.c

clean:
	rm -rf *.o GrishaginTest
