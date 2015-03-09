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

//Count the number of neighbors that has orientation q
//Used by rotChange
int Cell::numOfNNOr(int q) {
  int c = 0;

  if (im->orient == q) c++;
  if (ip->orient == q) c++;
  if (jm->orient == q) c++;
  if (jp->orient == q) c++;
  if (km->orient == q) c++;
  if (kp->orient == q) c++;

  return c;
}

//Count the number of neighbors that has id i
//Used by swapChange
int Cell::numOfNNId(int i) {
  int c = 0;

  if (im->id == i) c++;
  if (ip->id == i) c++;
  if (jm->id == i) c++;
  if (jp->id == i) c++;
  if (km->id == i) c++;
  if (kp->id == i) c++;

  return c;
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
void Cell::swapIdOr(Cell* c) { //TODO: Implement theta changes
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


//1D periodic boundary condition
int Simulation::wrap1d(int coord, int dir, int step) {
  coord += step;

  switch(dir) {
    case 0:
      if (coord < 0) coord += Lx;
      if (coord >= Lx) coord -= Lx;
      break;
    case 1:
      if (coord < 0) coord += Ly;
      if (coord >= Ly) coord -= Ly;
      break;
    case 2:
      if (coord < 0) coord += Lz;
      if (coord >= Lz) coord -= Lz;
      break;
  }

  return coord;
}
//Periodic boundary checking function
int Simulation::wrap3d(int index, int dir, int step) {
  //Convert index to 3D coords
  int k = index/(Lx*Ly); //Which slice
  int j = (index%(Lx*Ly))/Lx; //Which row on the slice
  int i = (index%(Lx*Ly))%Lx; //Which element in the row on the slice

  //Direction: 0 is i, 1 is j, 2 is k
  switch (dir) {
    case 0:
      wrap1d(i, dir, step);
      break;
    case 1:
      wrap1d(j, dir, step);
      break;
    case 2:
      wrap1d(k, dir, step);
      break;
  }

  //Revert back to index
  return i + Lx*j + Lx*Ly*k;
}

//Calculate the energy change for a rotation move
double Simulation::rotChange(Cell* c, int q) {
  //Look at how many orientations match in old and new config
  int o = c->numOfNNOr();
  int n = c->numOfNNOr(q);

  //cout << "o " << o << " n " << n << endl;

  //de = -A*(n - o)
  return (double)(o - n)*A;
}

//Calculate the energy change for a swap move
double Simulation::swapChange(Cell* c1, Cell* c2) {
  //Calculate current id and orientation matches
  int o_id1 = c1->numOfNNId(), o_id2 = c2->numOfNNId();
  int o_or1 = c1->numOfNNOr(), o_or2 = c2->numOfNNOr();


  //Get id and orientation info
  int i1 = c1->getId(), i2 = c2->getId();
  int o1 = c1->getOr(), o2 = c2->getOr();

  //Calculate proposed id and orientation matches
  int n_id1 = c1->numOfNNId(i2), n_id2 = c2->numOfNNId(i1);
  int n_or1 = c1->numOfNNOr(o2), n_or2 = c2->numOfNNOr(o1);

  //de = -K((n_id1 + n_id2) - (o_id1 + o_id2)) - A((n_or1 + n_or2) - (o_or1 + o_or2))
  double de = K*(double)((o_id1 + o_id2) - (n_id1 + n_id2)) + A*(double)((o_or1 + o_or2) - (n_or1 + n_or2));

  //If the two are neighbors and are diff id or orient, need to subtract an overcount
  //This is because the environment was not actually updated
  if (c1->isNeighbor(c2)) {
    if (i1 != i2) {
      de += 2*K;
    }
    if (o1 != o2) {
      de += 2*A;
    }
  }

  return de;
}

//Calculate theta (M=26)
//From V. Argawal and B. Peters, J. Chem. Phys. 140, 084111
double* Simulation::calctheta() {
  //Initialize
  double* theta = new double [NMAX];

  //Run through and calculate for every lattice position
  for (int i = 0; i < NMAX; i++) {
    theta[i] = 0.5;

    //Split into 3D coords
    int x = (i%(Lx*Ly))%Lx;
    int y = (i%(Lx*Ly))/Lx;
    int z = i/(Lx*Ly);

    //Run through neighbors (3x3x3 cube)
    int id = array[i]->getId();
    for (int a = -1; a < 2; a++) {
      for (int b = -1; b < 2; b++) {
        for (int c = -1; c < 2; c++) {
          //Convert back to 1D coord
          int index = wrap1d(x, 0, a) + wrap1d(y, 1, b)*Lx + wrap1d(z, 2, c)*Lx*Ly;

          //Increase theta if id match
          if (id == array[index]->getId()) {
            theta[i] += 1.0/52.0;
          }
        }
      }
    }

    //Remove overcount
    theta[i] -= 1.0/52.0;

    //If id = 2, we really wanted + to be -
    if (id == 2) {
      theta[i] = 1.0 - theta[i];
    }
  }

  return theta;
}

//Calculate Theta (M=27)
//From V. Argawal and B. Peters, J. Chem. Phys. 140, 084111
double* Simulation::calcTheta() {
  //Initialize
  double* Theta = new double [NMAX];

  //Get theta
  double* theta = calctheta();

  for (int i = 0; i < NMAX; i++) {
    Theta[i] = 0.0;

    //Split into 3D coords
    int x = (i%(Lx*Ly))%Lx;
    int y = (i%(Lx*Ly))/Lx;
    int z = i/(Lx*Ly);

    //Run through neighbors
    for (int a = -1; a <= 1; a++) {
      for (int b = -1; b <= 1; b++) {
        for (int c = -1; c <= 1; c++) {
          Theta[i] += theta[wrap1d(x, 0, a) + wrap1d(y, 1, b)*Lx + wrap1d(z, 2, c)*Lx*Ly];
        }
      }
    }

    //Average out Theta
    Theta[i] /= 27.0;
  }

  //Memory Management
  delete[] theta;

  return Theta;
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
    array[i]->setIm(array[wrap3d(i, 0, -1)]);
    array[i]->setIp(array[wrap3d(i, 0, 1)]);
    array[i]->setJm(array[wrap3d(i, 1, -1)]);
    array[i]->setJp(array[wrap3d(i, 1, 1)]);
    array[i]->setKm(array[wrap3d(i, 2, -1)]);
    array[i]->setKp(array[wrap3d(i, 2, 1)]);
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

  double* Theta = calcTheta();

  for (int i = 0; i < NMAX; i++) {
    //Get index in histogram
    int index = floorf(Theta[i]*100);
    if (index < 0) index = 0;
    if (index > 99) index = 99;

    h[index] += 1;
  }

  //Normalize
  for (int i = 0; i < 100; i++) {
    h[i] /= (double)NMAX;
  }

  //Memory Management
  delete[] Theta;

  return h;
}

double* Simulation::calcX1() {
  int n [2] = {0, 0}; //Number of particles in each phase
  int n1 [2] = {0, 0}; // Number of species i in each phase

  double* Theta = calcTheta();

  for (int i = 0; i < NMAX; i++) {
    //Decide if in 1-rich or 2-rich phase
    if (Theta[i] >= cutoff) {
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

  //Memory Management
  delete[] Theta;

  //Calculate mole fraction from n1/n for each phase
  double* X1 = new double [2];
  for (int i = 0; i < 2; i++) {
    X1[i] = (double)n1[i]/(double)n[i];
  }

  return X1;
}
