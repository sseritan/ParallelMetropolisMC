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
#include <sstream>
//Array include
#include <array>
//Stdlib include for srand
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
  double kT, compA, paramCutoff;
  if (argc < 6) {
    cout << "Usage:\n./MetropolisMC <Run Type> <Temp> <Eq. Sweeps> <Data Sweeps> <\%A> [<Cutoff>]" << endl;
    cout << "Run Types:" << endl;
    cout << "1) Melting Temperature. Only energy is kept." << endl;
    cout << "2) 2D Cutoff. Histogram of (Theta, Phi) is kept." << endl;
    cout << "3) Solid-Solid Phase Diagram. XA is calculated for two phases, using the cutoff value for Theta." << endl;
    cout << "4) Liquid-Solid Phase Diagram. XA is calculated for two phases, using the cutoff value for Phi." << endl;
    exit(1);
  } else {
    //Parse arguments
    runType = atoi(argv[1]);
    kT = (double)atof(argv[2]);
    eqSweeps = atoi(argv[3]);
    dataSweeps = atoi(argv[4]);
    compA = (double)atof(argv[5]);
    if (runType >= 3 && argc == 7) {
      paramCutoff = (double)atof(argv[6]);
    } else {
      cout << "No cutoff value specified." << endl;
      exit(1);
    }
  }

  //Print out simulation information
  cout << "SIMULATION DETAILS" << endl;
  cout << "Type of simulation: ";
  if (runType == 1) {
    cout << "Melting Temperature" << endl;
  } else if (runType == 2) {
    cout << "Cutoff (2D)" << endl;
  } else if (runType == 3) {
    cout << "Solid-Solid Phase Diagram" << endl;
  } else if (runType == 4) {
    cout << "Liquid-Solid Phase Diagram" << endl;
  } else {
    cout << "Invalid run type." << endl;
    exit(1);
  }

  auto start_time = chrono::system_clock::to_time_t(chrono::system_clock::now());
  cout << "Started simulation on " << ctime(&start_time) << endl;
  cout << "Box Dimensions: " << Lx << "x" << Ly << "x" << Lz << endl;
  cout << "Temperature: " << kT << endl;
  cout << "Hamiltonian: K=" << K << " A=" << A << endl;
  cout << "Move probabilities: Rotation=" << ROTATION << " Particle Swap=" << PARTSWAP << endl;
  cout << "Sweep information: Equilibration=" << eqSweeps << " Data Gathering=" << dataSweeps << endl;
  cout << "Seeding box with " << (int)floorf(compA*Lx*Ly*Lz) << " particles of species A" << endl;
  if (runType >= 3) {
    cout << "Using " << ((runType == 3) ? "Theta" : "Phi") << " cutoff value of " << paramCutoff << endl;
  }

  //Initialize random number generator
  auto seed = chrono::high_resolution_clock::now().time_since_epoch().count();
  cout << "Using seed " << seed << " for srand()." << endl;
  srand(seed);

  //Initialize identity and orientation arrays
  SimArray<int>* X_ptr = new SimArray<int>;
  SimArray<int>* S_ptr = new SimArray<int>;

  for(int i=0; i<Lx; i++) {
    for(int j=0; j<Ly; j++) {
      for(int k=0; k<Lz; k++) {
        //Seed identity so that target composition is met
        (*X_ptr)[i][j][k] = ((k<compA*Lz) ? 1 : 2);
        //Uniform orientation
        (*S_ptr)[i][j][k] = 1;
      }
    }
  }

  //Get initial energy of the system
  double e = energy(*X_ptr, *S_ptr)/kT;
  cout << "\nInitial energy (E/kT): " << e << endl;

  //Run through equilibration sweeps
  for (int t=1; t <= eqSweeps; t++) {
    //Run one full sweep
    sweep(*X_ptr, *S_ptr, e, kT);

    //Print out equilibration data
    if (t%100 == 0) {
      cout << "Equilibration sweep " << t << "/" << eqSweeps;

      //Print energy (updated and fresh)
      cout << " E/kT(updated)=" << e <<" E/kT(direct)=" << energy(*X_ptr, *S_ptr)/kT << endl;
    }
  }

  //Initialize data collection variables
  //MeltingTemp variable
  double eAvg;
  //Cutoff run variable
  Dim2Array h;
  //Phase Diagram run variable
  array<float, 2> XA;

  //Run through data sweeps (collecting data every sweep)
  cout << "*********************\nBeginning data sweeps\n*********************" << endl;
  for (int t = 1; t <= dataSweeps; t++) {
    //Run one full sweep
    sweep(*X_ptr, *S_ptr, e, kT);

    if (runType == 1) {
      //Melting Temp run

      if (t == 1) {
        //Initialize energy
        eAvg = e;
      } else {
        eAvg = (eAvg*(t-1) + e)/(float)t;
      }
    } else if (runType == 2) {
      //Cutoff run

      //Calculate phase and orientation parameters
      SimArray<float> Theta = phase_parameter(*X_ptr);
      SimArray<float> Phi = orientation_parameter(*S_ptr);

      if (t == 1) {
        //Get 2D histogram
        auto h = histogram2d(Theta, Phi);
      } else {
        //Update 2D histogram
        auto hNew = histogram2d(Theta, Phi);
        for (int i = 0; i < 100; i++) {
          for (int j = 0; j < 100; j++) {
            h[i][j] = (h[i][j]*(t-1) + hNew[i][j])/(float)t;
          }
        }
      }
    } else if (runType == 3 || runType == 4) {
      //Solid-Solid or Liquid-Solid Phase diagram run

      //Calculate phase parameters (based on type of run)
      SimArray<float> param = ((runType == 3) ? phase_parameter(*X_ptr) : orientation_parameter(*S_ptr));

      //Get phase composition data
      array<float, 2> XANew = phase_data(*X_ptr, param, paramCutoff);

      if (t == 1) {
        XA = XANew;
      } else {
        for (int i = 0; i < 2; i++) {
          XA[i] = (XA[i]*(t-1) + XANew[i])/(float)t;
        }
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
  if (runType == 1) {
    cout << "Average energy (E/kT): " << eAvg << endl;
  } else if (runType == 2) {
    cout << "2D Histogram:";
    for (auto a : h) {
      for (float f : a) {
        cout << f << " ";
      }
    }
    cout << endl;
  } else if (runType == 3 || runType == 4) {
    cout << "XA (A-rich): " << XA[0] << " XA (B-rich):" << XA[1] << endl;
  }

  auto end_time = chrono::system_clock::to_time_t(chrono::system_clock::now());
  cout << "\nFinished simulation on " << ctime(&end_time) << endl;

  //Memory Management
  delete X_ptr; delete S_ptr;

  //Exit successfully
  return 0;
}
