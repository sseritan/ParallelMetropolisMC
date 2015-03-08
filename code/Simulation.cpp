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

//Initialize theta (M=6) for a lattice position
//From V. Argawal and B. Peters, J. Chem. Phys. 140, 084111
void Cell::thetaInit() {
  theta = 0.5;

  //Run through neighbors
  if (id == im->id) {
    theta += 1.0/12.0;
  }
  if (id == ip->id) {
    theta += 1.0/12.0;
  }
  if (id == jm->id) {
    theta += 1.0/12.0;
  }
  if (id == jp->id) {
    theta += 1.0/12.0;
  }
  if (id == km->id) {
    theta += 1.0/12.0;
  }
  if (id == kp->id) {
    theta += 1.0/12.0;
  }

  //If id = 2, we actually wanted to flip it towards 0
  if (id == 2) {
    theta = 1.0 - theta;
  }
}

//Average thetas into Theta (M=7) for a lattice position
//From V. Argawal and B. Peters, J. Chem. Phys. 140, 084111
double Cell::calcTheta() {
  //Sum thetas
  double Theta = theta;
  Theta += im->theta;
  Theta += ip->theta;
  Theta += jm->theta;
  Theta += jp->theta;
  Theta += km->theta;
  Theta += kp->theta;

  //Average out Theta
  return Theta/7.0;
}

//Checks to see if a Cell is a neighbor
int Cell::isNeighbor(Cell* c) {
  if (im == c || ip == c || jm == c || jp == c || km == c || km == c) {
    return 1;
  }

  return 0;
}

//Swaps the identity and orientation of two cells
//This is cheaper than changing pointers and rearranging neighbor connections
void Cell::swapIdOr(Cell* c) {
  //Save my info
  int tempId = id, tempOr = orient;

  //Pull c's info
  id = c->id; orient = c->orient;

  //Put my info into c
  c->id = tempId; c->orient = tempOr;

  return;
}

/******************************
 * Simulation PRIVATE FUNCTIONS
 ******************************/
//Periodic boundary checking function
int Simulation::mod(int n) {
  while (n < 0) {
    n += NMAX;
  }
  while (n >= NMAX) {
    n -= NMAX;
  }

  return n;
}

//Calculate the energy change for a rotation move
double Simulation::rotChange(Cell* c, int q) {
  //Calculate current energy
  double e1 = c->pointEnergy(-1, -1);

  //Calculate fake energy
  double e2 = c->pointEnergy(-1, q);

  return (e2 - e1);
}

//Calculate the energy change for a swap move
double Simulation::swapChange(Cell* c1, Cell* c2) {
  //Calculate current energy
  double e1 = c1->pointEnergy(-1, -1) + c2->pointEnergy(-1, -1);

  int i1 = c1->getId(), i2 = c2->getId();
  int o1 = c1->getOr(), o2 = c2->getOr();

  //Calculate the fake energy when accepted
  double e2 = c1->pointEnergy(i2, o2) + c2->pointEnergy(i1, o1);

  //If the two are neighbors and are diff id or orient, need to subtract an overcount
  if (c1->isNeighbor(c2)) {
    if (i1 != i2) {
      e2 += 2*K;
    }
    if (o1 != o2) {
      e2 += 2*A;
    }
  }

  return (e2 - e1);
}



/*****************************
 * Simulation PUBLIC FUNCTIONS
 *****************************/

//Constructor
Simulation::Simulation(int x, int y, int z, double T, double compA, double c) {
  Lx = x; Ly = y; Lz = z;
  NMAX = Lx*Ly*Lz;
  kT = T;
  cutoff = c;
  cout << "\nINPUT DETAILS" << endl;
  cout << "Lattice dimensions of " << Lx << "x" << Ly << "x" << Lz << endl;
  cout << "Temperature (kT) of " << kT << endl;

  //Allocate array of Cells
  array = new Cell* [NMAX];

  //Initialize cells
  int threshold = compA*NMAX;
  cout << "Simulation seeded with " << threshold << "/" << NMAX << " particles of species 1" << endl;
  cout << "Using Theta cutoff value of " << cutoff << endl;

  //Set species 1
  for (int i = 0; i < threshold; i++) {
    array[i] = new Cell(1, 1);
  }
  //Set species 2
  for (int i = threshold; i < NMAX; i++) {
    array[i] = new Cell(2, 1);
  }

  //Secondary cell initialization
  for (int i = 0; i < NMAX; i++) {
    //Link cell with neighbors
    array[i]->setIm(array[mod(i-1)]);
    array[i]->setIp(array[mod(i+1)]);
    array[i]->setJm(array[mod(i-Lx)]);
    array[i]->setJp(array[mod(i+Lx)]);
    array[i]->setKm(array[mod(i-Lx*Ly)]);
    array[i]->setKp(array[mod(i+Lx*Ly)]);

    //Initialize theta
    array[i]->thetaInit();
  }

  //Initialize random number generation
  auto seed = chrono::high_resolution_clock::now().time_since_epoch().count();
  cout << "Using seed " << seed << " for srand()" << endl;
  srand(seed);

  //Initialize energy of the system
  energy = 0.0;
  for (int i = 0; i < NMAX; i++) {
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
  for (int i = 0; i < NMAX; i++) {
    delete array[i];
  }
  delete[] array;
}

//Function to evolve simulation by one sweep
void Simulation::doSweep() {
  for (int t = 1; t <= NMAX; t++) {
    double m = (double)rand()/(double)RAND_MAX;
    if (m < ROTATION) {
      //Rotation move
      //Pick a random location
      int index = rand()%(NMAX);

      //Pick random orientation
      int q = rand()%6 + 1;

      //Calculate energy change associated with rotation
      double de = rotChange(array[index], q)/kT;

      //Check acceptance
      if ((double)rand()/(double)RAND_MAX < exp(-de)) {
        //Update orientation and history
        array[index]->setOr(q);
        array[index]->pushUpdate(t);

        //Update energy
        energy += de;
      }
    } else {
      //Particle swap move
      int index1 = rand()%(NMAX);
      int index2 = rand()%(NMAX);

      //Calculate energy change associated with swap
      double de = swapChange(array[index1], array[index2])/kT;

      //Check acceptance
      if ((double)rand()/(double)RAND_MAX < exp(-de)) {
        //Swap cells
        array[index1]->swapIdOr(array[index2]);

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

double* Simulation::calcThetaHistogram() {
  //Initialize histogram
  double* h = new double [100];
  for (int i = 0; i < 100; i++) {
    h[i] = 0.0;
  }

  for (int i = 0; i < NMAX; i++) {
    //Calculate Theta
    double Theta = array[i]->calcTheta();

    //Get index in histogram (if Theta = 1.00, still in bin 99)
    int index = (floorf(Theta*100) <= 99 ? : 99);

    h[index] += 1;
  }

  return h;
}

double* Simulation::calcX1() {
  int n [2] = {0, 0}; //Number of particles in each phase
  int n1 [2] = {0, 0}; // Number of species i in each phase

  for (int i = 0; i < NMAX; i++) {
    //Calculate Theta
    double Theta = array[i]->calcTheta();

    //Decide if in 1-rich or 2-rich phase
    if (Theta >= cutoff) {
      n[0]++;
      if (array[i]->getId() == 1) {
        n1[0]++;
      }
    } else {
      n[1]++;
      if (array[i]->getId() == 1) {
        n1[1]++;
      }
    }
  }

  //Calculate mole fraction from n1/n for each phase
  double* X1 = new double [2];
  for (int i = 0; i < 2; i++) {
    X1[i] = (double)n1[i]/(double)n[i];
  }

  return X1;
}
