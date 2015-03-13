CC=icc
TBBDIR=/opt/intel/composer_xe_2015.2.132/tbb/
CFLAGS=-O2 -march=native -std=c++11 -ltbb -I $(TBBDIR)/include

MC: MetropolisMC.cpp Simulation.cpp Simulation.hpp CellMove.hpp
	$(CC) $(CFLAGS) MetropolisMC.cpp Simulation.cpp -o MC

test: test.cpp Simulation.hpp Simulation.cpp SimTest.cpp SimTest.hpp CellMove.hpp
	$(CC) -std=c++11 test.cpp SimTest.cpp Simulation.cpp -o test

clean:
	rm -f MC test
