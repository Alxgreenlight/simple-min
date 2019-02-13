CXX=g++
CXXFLAGS=-c -std=c++11
CC=gcc
CCFLAGS=-c

GKLStest: example.o ../solver/solver.o rnd_gen.o gkls.o
	$(CXX) -o GKLStest example.o ../solver/solver.o rnd_gen.o gkls.o
	
omp: example.o ../solver_parallel/solver_omp.o rnd_gen.o gkls.o
	$(CXX) -o GKLSTestOMP example.o ../solver_parallel/solver_omp.o rnd_gen.o gkls.o -fopenmp

th: example.o ../solver_parallel/solver_threads.o rnd_gen.o gkls.o
	$(CXX) -o GKLSTestTh example.o ../solver_parallel/solver_threads.o rnd_gen.o gkls.o -pthread

example.o: example.cpp
	$(CXX) $(CXXFLAGS) example.cpp

../solver/solver.o: ../solver/solver.cpp
	$(CXX) $(CXXFLAGS) -o ../solver/solver.o ../solver/solver.cpp
	
../solver_parallel/solver_omp.o: ../solver_parallel/solver_omp.cpp
	$(CXX) $(CXXFLAGS) ../solver_parallel/solver_omp.cpp -o ../solver_parallel/solver_omp.o -fopenmp

../solver_parallel/solver_threads.o: ../solver_parallel/solver_threads.cpp
	$(CXX) $(CXXFLAGS) ../solver_parallel/solver_threads.cpp -o ../solver_parallel/solver_threads.o -lpthread

gkls.o: gkls.c
	$(CC) $(CCFLAGS) gkls.c

rnd_gen.o: rnd_gen.c
	$(CC) $(CCFLAGS) rnd_gen.c

clean:
	rm -rf *.o ../solver/*.o GKLStest