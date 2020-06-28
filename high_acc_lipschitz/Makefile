CXX=g++
CXXFLAGS=-c -std=c++11 -fopenmp -Wall
CC=gcc
CCFLAGS=-c
INCDIR=-I../
GKLSDIR=../GKLStest

high_acc_L: high_acc.o gkls.o rnd_gen.o
	$(CXX) -o high_acc_L.exe high_acc.o gkls.o rnd_gen.o -fopenmp

high_acc.o: high_acc.cpp
	$(CXX) $(CXXFLAGS) $(INCDIR) high_acc.cpp

gkls.o: $(GKLSDIR)/gkls.c
	$(CC) $(CCFLAGS) $(GKLSDIR)/gkls.c

rnd_gen.o: $(GKLSDIR)/rnd_gen.c
	$(CC) $(CCFLAGS) $(GKLSDIR)/rnd_gen.c

clean:
	rm -rf *.o high_acc_L

rebuild: clean high_acc_L
