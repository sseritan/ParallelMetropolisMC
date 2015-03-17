//
// Originally written by Stefan Seritan for the Peters Group at UCSB
// Rewritten by Stefan Seritan and Wei Dai on 03/07/15
//
// Originally a sequential Metropolis Monte Carlo simulation
// Parallelized for CS 140 Winter 2015 Final Project
//


//I/O includes
#include <iostream>
//Stdlib include for atoi, atof
#include <stdlib.h>
//Time include
#include <chrono>
#include <time.h>
#include <omp.h>

//Local include
#include "./Simulation.hpp"

using namespace std;

int main(int argc, char* argv[]) {
	MPI_Init(&argc,&argv);

	int myrank;
	int nprocs;
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Status status;

  double start_time = MPI_Wtime();

  int runType, x, y, z, eqSweeps, dataSweeps, size, NMAX;
  double kT, compA, cutoff;
  if (argc != 9) {
    cout << "Usage: ./MC <kT> <x> <y> <z> <Eq. Sweeps> <Data Sweeps> <CompA> <Cutoff>" << endl;
    exit(1);
  } else {
    //Parse arguments
    kT = (double)atof(argv[1]);
    x = atoi(argv[2]);
    y = atoi(argv[3]);
    z = atoi(argv[4]);
    eqSweeps = atoi(argv[5]);
    dataSweeps = atoi(argv[6]);
    compA = (double)atof(argv[7]);
    cutoff = (double)atof(argv[8]);
    size = dataSweeps / nprocs;
    NMAX = x * y * z;
  }

  double start, end, initDuration, eqDuration, commDuration, dataDuration;
  start = MPI_Wtime();

  Simulation* sim = new Simulation(x, y, z, kT, compA, cutoff);

  if (myrank == 0) { // Master
    //Set output precision
    cout.precision(2); cout << fixed;

    //Print out simulation information
    cout << "Started simulation" << endl;

    cout << "Sweep information: Equilibration=" << eqSweeps << " Data Gathering=" << dataSweeps << endl;

    end = MPI_Wtime();
    //Initialization timing
    initDuration = end - start;

    start = MPI_Wtime();

    //Run equilibration sweeps
    for (int t = 1; t <= eqSweeps; t++) {
      sim->doSweep();

      if (eqSweeps >= 10 && t%(eqSweeps/10) == 0) {
        cout << "Equilibration sweep " << t << "/" << eqSweeps;
        cout << " E/kT=" << sim->getEnergy() << endl;
      }
    }

    end = MPI_Wtime();
    eqDuration = end - start;

    //Calculate Theta histogram just to make sure cutoff value is ok
    double* histogram = sim->calcThetaHistogram();
    cout << "Theta Histogram:" << endl;
    for (int i = 0; i < 100; i++) {
      cout << i << ": " << histogram[i] << " ";
      if (i > 0 && i%9 == 0) cout << endl;
    }
    cout << endl;

    //Memory Management
    delete[] histogram;

    start = MPI_Wtime();
  }

  // Send equillibrium state
  MPI_Bcast(sim->array, NMAX, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&sim->energy, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  end = MPI_Wtime();
  commDuration = end - start;

  start = MPI_Wtime();

  //Data collection variables
  double eAvg;
  double X1 [2];

  for (int t = 1; t <= size; t++) {
    sim->doSweep();

    //Update Theta and calculate phase compositions from Theta and cutoff
    double* X1New = sim->calcX1();

    if (t == 1) {
      //Initialize data
      eAvg = sim->getEnergy();
      X1[0] = X1New[0]; X1[1] = X1New[1];
    } else {
      //Keep running averages
      eAvg = (eAvg*(t-1) + sim->getEnergy())/(double)t;
      X1[0] = (X1[0]*(t-1) + X1New[0])/(double)t;
      X1[1] = (X1[1]*(t-1) + X1New[1])/(double)t;
    }

    //Memory Management
    delete[] X1New;
  }

  double GeAvg;
  double GX1[2];

  MPI_Reduce(&eAvg, &GeAvg, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&X1, &GX1, 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  GeAvg /= nprocs;
  GX1[0] /= nprocs;
  GX1[1] /= nprocs;

  if (myrank == 0) {
    end = MPI_Wtime();
    dataDuration = end - start;

    //Print collected data
    cout << "Average energy (E/kT): " << GeAvg << endl;
    cout << "X1 (1-rich): " << GX1[0] << " X1 (2-rich) " << GX1[1] << endl;

    cout << "\nInitialization time (s): " << initDuration << endl;
    cout << "Communication time (s): " << commDuration << endl;
    cout << "Equilibrium time (s): " << eqDuration << endl;
    cout << "Data time (s): " << dataDuration << endl;

    double end_time = MPI_Wtime();
    cout << "\nSimulation Wtime (s): " << end_time - start_time <<  endl;
  }

  //Memory Management
  delete sim;

	MPI_Finalize();

  //Exit successfully
  return 0;
}
