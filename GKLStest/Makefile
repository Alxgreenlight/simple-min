CXX=g++
CXXFLAGS=-c -std=c++11
CC=gcc
CCFLAGS=-c

GKLStest: example.o ../solver/solver.o rnd_gen.o gkls.o
	$(CXX) -o GKLStest example.o ../solver/solver.o rnd_gen.o gkls.o

example.o: example.cpp
	$(CXX) $(CXXFLAGS) example.cpp

../solver/solver.o: ../solver/solver.cpp
	$(CXX) $(CXXFLAGS) -o ../solver/solver.o ../solver/solver.cpp

gkls.o: gkls.c
	$(CC) $(CCFLAGS) gkls.c

rnd_gen.o: rnd_gen.c
	$(CC) $(CCFLAGS) rnd_gen.c

clean:
rm -rf *.o ../solver/*.o GKLStest