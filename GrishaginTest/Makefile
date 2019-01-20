GrishaginTest: example.o ../solver/solver.o vagris.o
	g++ example.o ../solver/solver.o vagris.o -o Grishagintest

example.o: example.cpp
	g++ -c example.cpp -std=c++11

../solver/solver.o: ../solver/solver.cpp
	g++ -c ../solver/solver.cpp -std=c++11 -o ../solver/solver.o

vagris.o: vagris.c
	g++ -c vagris.c

clean:
	rm -rf *.o ../solver/*.o GrishaginTesttest