CC=g++
CFLAGS=-O2 -march=native -std=c++11

MC: MetropolisMC.cpp MetropolisMC.hpp Simulation.cpp Simulation.hpp
	$(CC) $(CFLAGS) MetropolisMC.cpp Simulation.cpp -o MC

test: test.cpp Simulation.hpp Simulation.cpp
	$(CC) $(CFLAGS) test.cpp Simulation.cpp -o test

clean:
	rm MC test
