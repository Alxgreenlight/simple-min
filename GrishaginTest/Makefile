CXX=g++
CXXFLAGS=-c -std=c++11 -Wall -fopenmp
CC=gcc
CCFLAGS=-c

GrishaginTest: example.o vagris.o
	$(CXX) -o GrishaginTest example.o vagris.o -fopenmp

Rv: Rexample.o vagris.o
	$(CXX) -o R_GrishaginTest Rexample.o vagris.o -fopenmp

example.o: example.cpp
	$(CXX) $(CXXFLAGS) example.cpp

Rexample.o: Rexample.cpp
	$(CXX) $(CXXFLAGS) Rexample.cpp

vagris.o: vagris.c
	$(CC) $(CCFLAGS) vagris.c

clean:
	rm -rf *.o GrishaginTest R_GrishaginTest
