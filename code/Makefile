CC=mpicc
CFLAGS=-O3 -std=c++11 -openmp

MC: MetropolisMC.cpp Simulation.cpp Simulation.hpp MoveLock.hpp
	$(CC) $(CFLAGS) MetropolisMC.cpp Simulation.cpp -o MC

test: test.cpp Simulation.hpp Simulation.cpp SimTest.cpp SimTest.hpp MoveLock.hpp
	$(CC) $(CFLAGS) test.cpp SimTest.cpp Simulation.cpp CellMove.hpp -o test

clean:
	rm -f MC test
