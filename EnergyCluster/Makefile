CXX=g++
CXXFLAGS=-c -std=c++11 -fopenmp
CC=gcc
CCFLAGS=-c

EClust: EnergyCluster.o
	$(CXX) -o EClust EnergyCluster.o -fopenmp

EnergyCluster.o: EnergyCluster.cpp
	$(CXX) $(CXXFLAGS) EnergyCluster.cpp

clean:
	rm -rf *.o EClust
