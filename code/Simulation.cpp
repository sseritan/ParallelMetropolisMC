//
// Written by Stefan Seritan on 03/06/15 for CS 140
// All member functions for the Cell and Simulation classes
//

#include <iostream>
#include <vector>
#include <chrono>
#include <cmath>

//Local include
#include "Simulation.hpp"
#include "MetropolisMC.hpp"

using namespace std;

/***********************
 * Cell PUBLIC FUNCTIONS
 ***********************/

//Constructor
Cell::Cell(int i, int o) {
  id = i;
  orient = o;
  updateHistory.push_back(0);
}

//Setters
void Cell::pushUpdate(int m) {
  updateHistory.push_back(m);
}

void Cell::resetUpdate() {
  updateHistory.clear();
  updateHistory.push_back(0);
}

//Compute energy of a given id and orientation with 1 neighbor
double Cell::pairEnergy(int i, int o, int q) {
  //If passed negative, use actual Cell values
  if (i < 0) i = id;
  if (o < 0) o = orient;

  Cell* n;

  //Get neighbor
  switch (q) {
    case 1: n = im; break;
    case 2: n = ip; break;
    case 3: n = jm; break;
    case 4: n = jp; break;
    case 5: n = km; break;
    case 6: n = kp; break;
    default:
      cout << "Bad q pass to pairEnergy" << endl;
      exit(1);
      break;
  }

  double e = 0.0;
  if (i == n->id) e -= K;
  if (o == n->orient) e -= A;

  return e;
}

//Compute energy with all neighbors
//Pass neg values to use actual Cell values
//Otherwise, can be used to calculate "fake" energies
double Cell::pointEnergy(int i, int o) {
  //If neg values passed in, use actual Cell values
  if (i < 0) i = id;
  if (o < 0) o = orient;

  double e = 0.0;
  for (int q = 1; q <= 6; q++) {
    e += pairEnergy(i, o, q);
  }

  return e;
}

/******************************
 * Simulation PRIVATE FUNCTIONS
 ******************************/
//Periodic boundary checking function
int Simulation::mod(int n, int k) {
  while (n < 0) {
    n += k;
  }
  while (n >= k) {
    n -= k;
  }

  return n;
}

//TODO
double Simulation::rotChange(int i, int j, int k, int q) {
  return 0.0;
}

//TODO
double Simulation::swapChange(int i, int j, int k, int ii, int jj, int kk) {
  return 0.0;
}

/*****************************
 * Simulation PUBLIC FUNCTIONS
 *****************************/

//Constructor
Simulation::Simulation(int x, int y, int z, double T, double compA, double c) {
  Lx = x; Ly = y; Lz = z;
  kT = T;
  cutoff = c;
  cout << "\nINPUT DETAILS" << endl;
  cout << "Lattice dimensions of " << Lx << "x" << Ly << "x" << Lz << endl;
  cout << "Temperature (kT) of " << kT << endl;

  //Allocate array of Cells
  int nmax = Lx*Ly*Lz;
  array = new Cell* [nmax];

  //Initialize cells
  int threshold = compA*nmax;
  cout << "Simulation seeded with " << threshold << "/" << nmax << " particles of species 1" << endl;
  cout << "Using Theta cutoff value of " << cutoff << endl;

  //Set species 1
  for (int i = 0; i < threshold; i++) {
    array[i] = new Cell(1, 1);
  }
  //Set species 2
  for (int i = threshold; i < Lx*Ly*Lz; i++) {
    array[i] = new Cell(2, 1);
  }

  //Link array
  for (int i = 0; i < nmax; i++) {
    array[i]->setIm(array[mod(i-1, nmax)]);
    array[i]->setIp(array[mod(i+1, nmax)]);
    array[i]->setJm(array[mod(i-Lx, nmax)]);
    array[i]->setJp(array[mod(i+Lx, nmax)]);
    array[i]->setKm(array[mod(i-Lx*Ly, nmax)]);
    array[i]->setKp(array[mod(i+Lx*Ly, nmax)]);
  }

  //Initialize random number generation
  auto seed = chrono::high_resolution_clock::now().time_since_epoch().count();
  cout << "Using seed " << seed << " for srand()" << endl;
  srand(seed);

  //Initialize energy of the system
  energy = 0.0;
  for (int i = 0; i < nmax; i++) {
    //Calculate pairwise energy with forward pairs in each direction
    energy += array[i]->pairEnergy(-1, -1, 2); //In forward i direction
    energy += array[i]->pairEnergy(-1, -1, 4); //In forward j direction
    energy += array[i]->pairEnergy(-1, -1, 6); //In forward k direction
  }

  cout << "Initial energy: " << energy/kT << endl;
  cout << endl;
}

//Destructor
Simulation::~Simulation() {
  //Memory Management
  for (int i = 0; i < Lx*Ly*Lz; i++) {
    delete array[i];
  }
  delete[] array;
}

//Function to evolve simulation by one sweep
void Simulation::doSweep() {
  for (int t = 1; t <= Lx*Ly*Lz; t++) {
    double m = (double)rand()/(double)RAND_MAX;
    if (m < ROTATION) {
      //Rotation move
      int i = rand()%Lx;
      int j = rand()%Ly;
      int k = rand()%Lz;

      //Pick random orientation
      int q = rand()%6 + 1;

      //Calculate energy change associated with rotation
      double de = rotChange(i, j, k, q)/kT;

      if ((double)rand()/(double)RAND_MAX < exp(-de)) {
        //Accept move
        int index = i + Lx*j + Lx*Ly*k;

        //Update orientation and history
        array[index]->setOr(q);
        array[index]->pushUpdate(t);

        //Update energy
        energy += de;
      }
    } else {
      //Particle swap move
      int i = rand()%Lx; int ii = rand()%Lx;
      int j = rand()%Ly; int jj = rand()%Ly;
      int k = rand()%Lz; int kk = rand()%Lz;

      //Calculate energy change associated with swap
      double de = swapChange(i, j, k, ii, jj, kk)/kT;

      //Check acceptance
      if ((double)rand()/(double)RAND_MAX < exp(-de)) {
        int index1 = i + Lx*j + Lx*Ly*k;
        int index2 = ii + Lx*jj + Lx*Ly*kk;

        //Swap cell pointers
        Cell* temp = array[index1];
        array[index1] = array[index2];
        array[index2] = temp;

        //Update history
        array[index1]->pushUpdate(t);
        array[index2]->pushUpdate(t);

        //Update energy
        energy += de;
      }
    }
  }
}

//Function to return energy
double Simulation::getEnergy() {
  return energy;
}

void Simulation::updateTheta() {
  //TODO
}

double* Simulation::calcThetaHistogram() {
  //TODO
  return NULL;
}

double* Simulation::calcX1() {
  //TODO
  return NULL;
}
