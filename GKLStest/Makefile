GKLStest: example.o ../solver/solver.o rnd_gen.o gkls.o
	g++ example.o ../solver/solver.o rnd_gen.o gkls.o -lm -o GKLStest

example.o: example.cpp
	g++ -c example.cpp -std=c++11

../solver/solver.o: ../solver/solver.cpp
	g++ -c ../solver/solver.cpp -std=c++11 -o ../solver/solver.o

gkls.o: gkls.c
	g++ -c gkls.c

rnd_gen.o: rnd_gen.c
	g++ -c rnd_gen.c

clean:
	rm -rf *.o ../solver/*.o GKLStest