CXX=g++
CXXFLAGS=-c -std=c++11
CC=gcc
CCFLAGS=-c

DiRectTest: tstc.o ../solver/solver.o rnd_gen.o gkls.o DIRect.o direct_wrap.o DIRserial.o DIRsubrout.o
$(CXX) -o DiRectTest tstc.o ../solver/solver.o rnd_gen.o gkls.o DIRect.o direct_wrap.o DIRserial.o DIRsubrout.o

example.o: tstc.cpp
	$(CXX) $(CXXFLAGS) example.cpp

../solver/solver.o: ../solver/solver.cpp
	$(CXX) $(CXXFLAGS) -o ../solver/solver.o ../solver/solver.cpp

gkls.o: gkls.c
	$(CC) $(CCFLAGS) gkls.c

rnd_gen.o: rnd_gen.c
	$(CC) $(CCFLAGS) rnd_gen.c

DIRect.o: DIRect.c
	$(CC) $(CCFLAGS) DIRect.c

direct_wrap.o: direct_wrap.c
	$(CC) $(CCFLAGS) direct_wrap.c

DIRserial.o: DIRserial.c
	$(CC) $(CCFLAGS) DIRserial.c

DIRsubrout.o: DIRsubrout.c
	$(CC) $(CCFLAGS) DIRsubrout.c

clean:
	rm -rf *.o ../solver/*.o DiRectTest