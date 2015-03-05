//
//
// Originally written by Stefan Seritan on 8/28/14
// Modified by Stefan Seritan and Wei Dai for CS 140 Final Project Winter 2015
//
// Originally a sequential Metropolis Monte Carlo simulation
// Parallelized for CS 140 Final Project
//


/***************************************************/
/*        HYBRID POTTS-LATTICE GAS MODEL           */
/***************************************************/

//I/O includes
#include <iostream>
//Array include
#include <array>
//Stdlib include for atoi, atof
#include <cstdlib>
//Math include
#include <cmath>
//Time include
#include <chrono>
#include <ctime>

//Local include
#include "./MetropolisMC.hpp"
#include "./subroutine.cpp"

using namespace std;

int main(int argc, char* argv[]) {
  //Set output precision
  cout.precision(2); cout << fixed;

  int runType, eqSweeps, dataSweeps;
  double kT, compA, cutoff;
  if (argc != 6) {
    cout << "Usage:\n./MetropolisMC <Temp> <Eq. Sweeps> <Data Sweeps> <CompA> <Cutoff>" << endl;
    cout << "Temp is kT. Acceptance is based on exp(-energy/kT). Higher kT basically allows less favored moves still get accepted, equivalent to higher temperature." << endl;
    cout << "Eq. Sweeps is the number of sweeps to reach equilibrium (no data collected). Equilibrium reached when energy doesn't change." << endl;
    cout << "Data. Sweeps is the number of sweeps where we collect energy and phase data between sweeps. More data sweeps, the better the collected statistics (weak law of large numbers)." << endl;
    cout << "CompA is a number between 0 and 1 that sets an initial fraction of the identities to 1 (species A)." << endl;
    cout << "Cutoff is the number we use to separate phases. Theta values above cutoff are one phase, values below are the other phase.\n"
     "For the symmetric Hamiltonian, it should be the same value as CompA." << endl;
    exit(1);
  } else {
    //Parse arguments
    kT = (double)atof(argv[1]);
    eqSweeps = atoi(argv[2]);
    dataSweeps = atoi(argv[3]);
    compA = (double)atof(argv[4]);
    cutoff = (double)atof(argv[5]);
  }

  //Print out simulation information
  cout << "SIMULATION DETAILS" << endl;
  auto start_time = chrono::system_clock::to_time_t(chrono::system_clock::now());
  cout << "Started simulation on " << ctime(&start_time) << endl;
  cout << "Box Dimensions: " << Lx << "x" << Ly << "x" << Lz << endl;
  cout << "Temperature: " << kT << endl;
  cout << "Hamiltonian: K=" << K << " A=" << A << endl;
  cout << "Move probabilities: Rotation=" << ROTATION << " Particle Swap=" << PARTSWAP << endl;
  cout << "Sweep information: Equilibration=" << eqSweeps << " Data Gathering=" << dataSweeps << endl;
  cout << "Seeding box with " << (int)floorf(compA*Lx*Ly*Lz) << " particles of species A" << endl;
  cout << "Using Theta cutoff value of " << cutoff << endl;

  //Initialize random number generator
  auto seed = chrono::high_resolution_clock::now().time_since_epoch().count();
  cout << "Using seed " << seed << " for srand()." << endl;
  srand(seed);

  //Initialize identity and orientation arrays
  SimArray<int>* X_ptr = new SimArray<int>(Lx, Ly, Lz);
  SimArray<int>* S_ptr = new SimArray<int>(Lx, Ly, Lz);

  //Initialize array for checking conflicts
  vector<Move>* M = new vector<Move>[np];

  for(int i=0; i<Lx; i++) {
    for(int j=0; j<Ly; j++) {
      for(int k=0; k<Lz; k++) {
        //Seed identity so that target composition is met
        (*X_ptr)(i,j,k) = ((k<compA*Lz) ? 1 : 2);
        //Uniform orientation
        (*S_ptr)(i,j,k) = 1;
      }
    }
  }

  //Get initial energy of the system
  double e = energy(*X_ptr, *S_ptr)/kT;
  cout << "\nInitial energy (E/kT): " << e << endl;

  //Run through equilibration sweeps
  for (int t=1; t <= eqSweeps; t++) {
    //Run one full sweep
    sweep(M, *X_ptr, *S_ptr, e, kT);

    //Print out equilibration data
    if (t%100 == 0) {
      cout << "Equilibration sweep " << t << "/" << eqSweeps;

      //Print energy (updated and fresh)
      cout << " E/kT(updated)=" << e << endl;
    }
  }

  //Initialize data collection variables
  //MeltingTemp variable
  double eAvg;
  //Phase Diagram run variable
  array<double, 2> XA;

  //Run through data sweeps (collecting data every sweep)
  cout << "*********************\nBeginning data sweeps\n*********************" << endl;
  for (int t = 1; t <= dataSweeps; t++) {
    //Run one full sweep
    sweep(M, *X_ptr, *S_ptr, e, kT);

    //Calculate Theta
    SimArray<double> Theta = phase_parameter(*X_ptr);

    //Calculate phase compositions from Theta and cutoff value
    array<double, 2> XANew = phase_data(*X_ptr, Theta, cutoff);

    if (t == 1) {
      //Initialize data
      eAvg = e;
      XA = XANew;
    } else {
      //Keep running avg of data
      eAvg = (eAvg*(t-1) + e)/(double)t;
      for (int i = 0; i < 2; i++) {
        XA[i] = (XA[i]*(t-1) + XANew[i])/(double)t;
      }
    }

    /*
    //Print snapshots every 10000 sweeps, just to see what's up
    if (t%10000 == 0) {
      cout << "Printing VMD snapshot at " << t << " data sweeps." << endl;
      print_VMD_snapshot(*X_ptr, *S_ptr, t);
    }
    */
  }

  //Print out information collected in data sweeps
  cout << "Average energy (E/kT): " << eAvg << endl;
  cout << "XA (A-rich): " << XA[0] << " XA (B-rich):" << XA[1] << endl;

  auto end_time = chrono::system_clock::to_time_t(chrono::system_clock::now());
  cout << "\nFinished simulation on " << ctime(&end_time) << endl;

  //Memory Management
  delete X_ptr; delete S_ptr; delete[] M;

  //Exit successfully
  return 0;
}
