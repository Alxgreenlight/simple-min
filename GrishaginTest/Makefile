CXX=g++
CXXFLAGS=-c -std=c++11 -Wall
CC=gcc
CCFLAGS=-c

GrishaginTest: example.o ../solver/solver.o vagris.o
	$(CXX) -o GrTest example.o ../solver/solver.o vagris.o

example.o: example.cpp
	$(CXX) $(CXXFLAGS) example.cpp

../solver/solver.o: ../solver/solver.cpp
	$(CXX) $(CXXFLAGS) ../solver/solver.cpp -o ../solver/solver.o

vagris.o: vagris.c
	$(CC) $(CCFLAGS) vagris.c

clean:
	rm -rf *.o ../solver/*.o GrTest