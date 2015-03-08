CC=g++
CFLAGS=-O2 -march=native -std=c++11

MC: MetropolisMC.cpp MetropolisMC.hpp Simulation.cpp Simulation.hpp
	$(CC) $(CFLAGS) MetropolisMC.cpp Simulation.cpp -o MC

MetropolisMC: MetropolisMC.cpp MetropolisMC.hpp subroutine.cpp
	$(CC) $(CFLAGS) MetropolisMC.cpp -o MetropolisMC

clean:
	rm MC
