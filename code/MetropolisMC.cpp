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
#include <cstdlib>
//Time include
#include <chrono>
#include <ctime>

//Local include
#include "./MetropolisMC.hpp"
#include "./Simulation.hpp"

using namespace std;

int main(int argc, char* argv[]) {
  //Set output precision
  cout.precision(2); cout << fixed;

  int runType, x, y, z, eqSweeps, dataSweeps;
  double kT, compA, cutoff;
  if (argc != 9) {
    cout << "Usage:\n./MC <kT> <x> <y> <z> <Eq. Sweeps> <Data Sweeps> <CompA> <Cutoff>" << endl;
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
  }

  //Print out simulation information
  auto start_time = chrono::system_clock::to_time_t(chrono::system_clock::now());
  cout << "Started simulation on " << ctime(&start_time) << endl;
  cout << "SIMULATION DETAILS" << endl;
  cout << "Hamiltonian: K=" << K << " A=" << A << endl;
  cout << "Move probabilities: Rotation=" << ROTATION << " Particle Swap=" << PARTSWAP << endl;
  cout << "Sweep information: Equilibration=" << eqSweeps << " Data Gathering=" << dataSweeps << endl;

  Simulation* sim = new Simulation(x, y, z, kT, compA, cutoff);

  //Run equilibration sweeps
  for (int t = 1; t <= eqSweeps; t++) {
    sim->doSweep();

    if (eqSweeps >= 10 && t%(eqSweeps/10) == 0) {
      cout << "Equilibration sweep " << t << "/" << eqSweeps;
      cout << " E/kT=" << sim->getEnergy() << endl;
    }
  }

  //Calculate Theta histogram just to make sure cutoff value is ok
  sim->updateTheta();
  double* histogram = sim->calcThetaHistogram();
  cout << "Theta Histogram:" << endl;
  for (int i = 0; i < 100; i++) {
    cout << "0: " << histogram[i] << " ";
    if (i%10 == 0) cout << endl;
  }
  cout << endl;

  //Data collection variables
  double eAvg;
  double X1 [2];

  for (int t = 1; t <= dataSweeps; t++) {
    sim->doSweep();

    //Update Theta and calculate phase compositions from Theta and cutoff
    sim->updateTheta();
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

    if (dataSweeps >= 10 && t%(dataSweeps/10) == 0) {
      cout << "Data sweep " << t << "/" << dataSweeps << endl;
    }
  }

  //Print collected data
  cout << "Average energy (E/kT): " << eAvg << endl;
  cout << "X1 (1-rich):" << X1[0] << " X1 (2-rich) " << X1[1] << endl;

  auto end_time = chrono::system_clock::to_time_t(chrono::system_clock::now());
  cout << "\nFinished simulation on " << ctime(&end_time) << endl;

  //Memory Management
  delete sim;

  //Exit successfully
  return 0;
}
