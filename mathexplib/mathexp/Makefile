CXX=g++
CXXFLAGS=-c -std=c++11

mathexp: main.o ../../solver/solver.o
$(CXX) -o mathexp main.o ../../solver/solver.o

omp: main.o ../../solver_parallel/solver_omp.o
$(CXX) -o mathexpOMP main.o ../../solver_parallel/solver_omp.o -fopenmp

main.o: main.cpp
$(CXX) $(CXXFLAGS) main.cpp

../../solver/solver.o: ../../solver/solver.cpp
$(CXX) $(CXXFLAGS) -o ../../solver/solver.o ../../solver/solver.cpp

../../solver_parallel/solver_omp.o: ../../solver_parallel/solver_omp.cpp
$(CXX) $(CXXFLAGS) ../../solver_parallel/solver_omp.cpp -o ../../solver_parallel/solver_omp.o -fopenmp

clean:
rm -rf *.o mathexp mathexpOMP