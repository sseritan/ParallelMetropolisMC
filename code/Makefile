CC=g++-4.9
CFLAGS=-O3 -std=c++11 -I /usr/local/include -L /usr/local/lib -ltbb

MC: MetropolisMC.cpp Simulation.cpp Simulation.hpp CellMove.hpp
	$(CC) $(CFLAGS) MetropolisMC.cpp Simulation.cpp CellMove.hpp -o MC

test: test.cpp Simulation.hpp Simulation.cpp SimTest.cpp SimTest.hpp CellMove.hpp
	$(CC) $(CFLAGS) test.cpp SimTest.cpp Simulation.cpp CellMove.hpp -o test

clean:
	rm -f MC test
