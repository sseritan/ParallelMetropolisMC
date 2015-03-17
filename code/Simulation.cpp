//
// Written by Stefan Seritan on 03/06/15 for CS 140
// All member functions for the Cell and Simulation classes
//

#include <iostream>
#include <chrono>
#include <cmath>

//Local include
#include "./Simulation.hpp"

using namespace std;

/******************************
 * Simulation PRIVATE FUNCTIONS
 ******************************/


//Calculate the energy change for a rotation move
double Simulation::rotChange(int pos, int q) {
  //Look at how many orientations match in old and new config
  int o = numOfNNOr(pos, array[pos]%10);
  int n = numOfNNOr(pos, q);

  //de = -A*(n - o)
  return (double)(o - n)*A;
}

//Calculate the energy change for a swap move
double Simulation::swapChange(int pos1, int pos2) {
  //Get id and orientation info
  int c1 = array[pos1], c2 = array[pos2];

  //Calculate difference in id matches
  int did = (numOfNNId(pos1, c2/10) + numOfNNId(pos2, c1/10)) - (numOfNNId(pos1, c1/10) + numOfNNId(pos2, c2/10));

  //Calculate difference in orientation matches
  int dor = (numOfNNOr(pos1, c2%10) + numOfNNOr(pos2, c1%10)) - (numOfNNOr(pos1, c1%10) + numOfNNOr(pos2, c2%10));

  //de = -K((n_id1 + n_id2) - (o_id1 + o_id2)) - A((n_or1 + n_or2) - (o_or1 + o_or2))
  double de = -K*(double)did - A*(double)dor;

  //If the two are neighbors and are diff id or orient, need to subtract an overcount
  //This is because the environment was not actually updated
  if (areNN(pos1, pos2)) {
    if (c1/10 != c2/10) {
      de += 2*K;
    }
    if (c1%10 != c2%10) {
      de += 2*A;
    }
  }

  return de;
}

void Simulation::performMove(Move* m) {
  int type = m->type;
  int pos = m->pos;
  int par = m->par;

  if (!type) {
    //Get energy change associated with rotation
    double de = rotChange(pos, par)/kT;

    //Check acceptance
    if ((double)rand()/(double)RAND_MAX < exp(-de)) {
      //Update orientation and history
      int curr = array[pos]%10;
      array[pos] += par - curr;

      //Update energy
      energy += de;
    }
  } else {
    //Calculate energy change associated with swap
    double de = swapChange(pos, par)/kT;

    //Check acceptance
    if ((double)rand()/(double)RAND_MAX < exp(-de)) {
      //Swap cells
      int temp = array[pos];
      array[pos] = array[par];
      array[par] = temp;

      //Update energy
      energy += de;
    }
  }
}

//Calculate theta (M=26)
//From V. Argawal and B. Peters, J. Chem. Phys. 140, 084111
double* Simulation::calctheta() {
  //Initialize
  double* theta = new double [NMAX];

  //Run through and calculate for every lattice position
  //TODO: parallel_for (blocks?)
  for (int i = 0; i < NMAX; i++) {
    theta[i] = 0.5;

    //Split into 3D coords
    int x = (i%(Lx*Ly))%Lx;
    int y = (i%(Lx*Ly))/Lx;
    int z = i/(Lx*Ly);

    //Run through neighbors (3x3x3 cube)
    int id = array[i]/10;
    for (int a = -1; a < 2; a++) {
      for (int b = -1; b < 2; b++) {
        for (int c = -1; c < 2; c++) {
          //Convert back to 1D coord
          int index = wrap1d(x, 0, a) + wrap1d(y, 1, b)*Lx + wrap1d(z, 2, c)*Lx*Ly;

          //Increase theta if id match
          if (id == array[index]/10) {
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

  //TODO: parallel_for (blocks?)
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

//Calculate the energy between two lattice positions
double Simulation::pairEnergy(int pos1, int pos2) {
  double e = 0.0;

  if (array[pos1]/10 == array[pos2]/10) e -= K;
  if (array[pos1]%10 == array[pos2]%10) e -= A;

  return e;
}

//Check to see how many neighbors of array[pos] have id i
int Simulation::numOfNNId(int pos, int i) {
  int count = 0;

  //Run through neighbors and count matching ids
  for (int j = 0; j < 3; j++) {
    for (int k = -1; k <= 1; k += 2) {
      if (array[step3d(pos, j, k)]/10 == i) count++;
    }
  }

  return count;
}

//Check to see how many neighbors of array[pos] have orientation o
int Simulation::numOfNNOr(int pos, int o) {
  int count = 0;

  //Run through neighbors and count matching ids
  for (int i = 0; i < 3; i++) {
    for (int j = -1; j <= 1; j += 2) {
      if (array[step3d(pos, i, j)]%10 == o) count++;
    }
  }

  return count;
}

//Check to see if two indices are neighbors in 3D
int Simulation::areNN(int pos1, int pos2) {
  for (int i = 0; i < 3; i++) {
    for (int j = -1; j <= 1; j+= 2) {
      if (step3d(pos1, i, j) == pos2) return 1;
    }
  }

  return 0;
}

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

//Step of 1D index in 3D coords, with periodic boundary conditions
int Simulation::step3d(int index, int dir, int step) {
  //Convert index to 3D coords
  int k = index/(Lx*Ly); //Which slice
  int j = (index%(Lx*Ly))/Lx; //Which row on the slice
  int i = (index%(Lx*Ly))%Lx; //Which element in the row on the slice

  //Direction: 0 is i, 1 is j, 2 is k
  switch (dir) {
    case 0:
      i = wrap1d(i, dir, step);
      break;
    case 1:
      j = wrap1d(j, dir, step);
      break;
    case 2:
      k = wrap1d(k, dir, step);
      break;
  }

  //Revert back to index
  return i + Lx*j + Lx*Ly*k;
}


/*****************************
 * Simulation PUBLIC FUNCTIONS
 *****************************/


//Constructor
Simulation::Simulation(int x, int y, int z, double T, double compA, double c) : Lx(x), Ly(y), Lz(z), NMAX(x*y*z), kT(T), cutoff(c) {
  cout << "Hamiltonian: K=" << K << " A=" << A << endl;
  cout << "Move probabilities: Rotation=" << ROTATION << " Particle Swap=" << PARTSWAP << endl;
  cout << "Lattice dimensions of " << Lx << "x" << Ly << "x" << Lz << endl;
  cout << "Temperature (kT) of " << kT << endl;

  //Allocate array (10's place is id, 1's place is orientation)
  array = new int [NMAX];

  //Initialize cells
  int threshold = compA*NMAX;
  cout << "Simulation seeded with " << threshold << "/" << NMAX << " particles of species 1" << endl;
  cout << "Using Theta cutoff value of " << cutoff << endl;

  //Set species 1
  for (int i = 0; i < threshold; i++) {
    array[i] = 11;
  }
  //Set species 2
  for (int i = threshold; i < NMAX; i++) {
    array[i] = 21;
  }

  //Initialize random number generation
  auto seed = chrono::high_resolution_clock::now().time_since_epoch().count();
  cout << "Using seed " << seed << " for srand()" << endl;
  srand(seed);

  //Initialize energy of the system
  energy = 0.0;
  for (int i = 0; i < NMAX; i++) {
    //Calculate pairwise energy with forward pairs in each direction
    for (int j = 0; j < 3; j++) {
      energy += pairEnergy(i, step3d(i, j, 1));
    }
  }
  energy /= kT;

  cout << "Initial energy: " << energy << endl;
  cout << endl;
}

//Destructor
Simulation::~Simulation() {
  //Memory Management
  delete[] array;
}

//Function to evolve simulation by one sweep
void Simulation::doSweep() {
  //Generate moves into array
  Move* moves [NMAX];
  for (int i = 0; i < NMAX; i++) {
    //Decide rotation (0) or swap (1)
    int type = (((double)rand()/(double)RAND_MAX < ROTATION) ? 0 : 1);
    int pos = rand()%NMAX;
    int param = (type ? rand()%NMAX : rand()%6 + 1); //New orient if rot, 2nd lattice position if swap
    moves[i] = new Move(type, pos, param);
  }

  //Run through moves and perform them (in serial bcuz overhead not worth it)
  for (int i = 0; i < NMAX; i++) {
    performMove(moves[i]);

    //Memory Managament
    delete moves[i];
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
      if (array[i]/10 == 1) {
        n1[0]++;
      }
    } else {
      n[1]++;
      if (array[i]/10 == 1) {
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
