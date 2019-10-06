CXX=g++
CXXFLAGS=-c -std=c++11 -fopenmp -Wall

miniOneDim: OneDimStarter.o
	$(CXX) -o miniOneDim OneDimStarter.o -fopenmp

OneDimStarter.o: OneDimStarter.cpp
	$(CXX) $(CXXFLAGS) OneDimStarter.cpp

clean:
	rm -rf *.o miniOneDim
