CXX=g++
CXXFLAGS=-c -std=c++11 -O -fopenmp -Wall
INCDIR=-I../

miniOneDim: OneDimStarter.o
	$(CXX) -o miniOneDim -O OneDimStarter.o -fopenmp

OneDimStarter.o: OneDimStarter.cpp
	$(CXX) $(CXXFLAGS) $(INCDIR) OneDimStarter.cpp

clean:
	rm -rf *.o miniOneDim

rebuild: clean miniOneDim
