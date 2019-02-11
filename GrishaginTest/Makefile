CXX=g++
CXXFLAGS=-c -std=c++11 -Wall
CC=gcc
CCFLAGS=-c

GrishaginTest: example.o ../solver/solver.o vagris.o
	$(CXX) -o GrTest example.o ../solver/solver.o vagris.o

omp: example.o ../solver_parallel/solver_omp.o vagris.o
	$(CXX) -o GrTestOMP example.o ../solver_parallel/solver_omp.o vagris.o -fopenmp

th: example.o ../solver_parallel/solver_threads.o vagris.o
	$(CXX) -o GrTestTh example.o ../solver_parallel/solver_threads.o vagris.o -pthread

example.o: example.cpp
	$(CXX) $(CXXFLAGS) example.cpp

../solver/solver.o: ../solver/solver.cpp
	$(CXX) $(CXXFLAGS) ../solver/solver.cpp -o ../solver/solver.o

../solver_parallel/solver_omp.o: ../solver_parallel/solver_omp.cpp
	$(CXX) $(CXXFLAGS) ../solver_parallel/solver_omp.cpp -o ../solver_parallel/solver_omp.o -fopenmp

../solver_parallel/solver_threads.o: ../solver_parallel/solver_threads.cpp
	$(CXX) $(CXXFLAGS) ../solver_parallel/solver_threads.cpp -o ../solver_parallel/solver_threads.o -lpthread

vagris.o: vagris.c
	$(CC) $(CCFLAGS) vagris.c

clean:
	rm -rf *.o ../solver/*.o ../solver_parallel/*.o GrTest GrTestOMP GrTestTh