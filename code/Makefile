CC=icc
CFLAGS=-O3 -std=c++11 -openmp

MetropolisMC: MetropolisMC.cpp MetropolisMC.hpp subroutine.cpp
	$(CC) $(CFLAGS) MetropolisMC.cpp -o MetropolisMC

clean:
	rm MetropolisMC
