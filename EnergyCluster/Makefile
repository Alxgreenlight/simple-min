CXX=g++
CXXFLAGS=-c -std=c++11 -O -fopenmp
CC=gcc
CCFLAGS=-c
INCDIR=-I../

EClust: EnergyCluster.o
	$(CXX) -o EClust -O EnergyCluster.o -fopenmp

EnergyCluster.o: EnergyCluster.cpp
	$(CXX) $(INCDIR) $(CXXFLAGS) EnergyCluster.cpp

clean:
	rm -rf *.o EClust

rebuild: clean EClust
