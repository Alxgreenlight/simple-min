CXX=g++
CXX=g++
CXXFLAGS=-c -std=c++11 -O -fopenmp
INCDIR=-I../../

testfuncs: testfuncs.o
	$(CXX) -o testfuncs testfuncs.o -O -fopenmp

testfuncs.o: testfuncs.cpp
	$(CXX) $(CXXFLAGS) $(INCDIR) testfuncs.cpp

clean:
	rm -rf *.o testfuncs

rebuild: clean testfuncs
