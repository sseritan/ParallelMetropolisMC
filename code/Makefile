CC=icc
CFLAGS=-O3 -march=native -std=c++11

MC: MetropolisMC.cpp Simulation.cpp Simulation.hpp CellMove.hpp
	$(CC) $(CFLAGS) MetropolisMC.cpp Simulation.cpp -o MC

test: test.cpp Simulation.hpp Simulation.cpp SimTest.cpp SimTest.hpp CellMove.hpp
	$(CC) -std=c++11 test.cpp SimTest.cpp Simulation.cpp -o test

clean:
	rm -f MC test
