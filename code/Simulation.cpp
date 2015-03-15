//
// Written by Stefan Seritan on 03/06/15 for CS 140
// Modified by Wei Dai 03/15/15
// All member functions for the Cell and Simulation classes
//

#include <iostream>
#include <vector>
#include <chrono>
#include <cmath>
#include <tbb/parallel_for.h>
#include <tbb/parallel_reduce.h>

//Local include
#include "./Simulation.hpp"
#include "./CellMove.hpp"

using namespace std;

/******************************
 * Simulation PRIVATE FUNCTIONS
 ******************************/

void Simulation::posLocks(int pos, int& even, int& odd) const {
  /* if (pos == 0) cout << "pos == 0" << endl; */
  //Split into 3D coords
  int x = (pos%(Lx*Ly))%Lx;
  int y = (pos%(Lx*Ly))/Lx;
  int z = pos/(Lx*Ly);

  /* cout << "x: " << x << endl; */
  /* cout << "y: " << y << endl; */
  /* cout << "z: " << z << endl; */

  /* cout << "even x: " << wrap1d(x, 0, -(x&1)) << endl; */
  /* cout << "even y: " << wrap1d(y, 1, -(y&1)) << endl; */
  /* cout << "even z: " << wrap1d(z, 2, -(z&1)) << endl; */

  /* cout << "odd x: " << wrap1d(x, 0, -((x&1)^1)) << endl; */
  /* cout << "odd y: " << wrap1d(y, 1, -((y&1)^1)) << endl; */
  /* cout << "odd z: " << wrap1d(z, 2, -((z&1)^1)) << endl; */

  //lock even lock
  even = wrap1d(x, 0, -(x&1)) + wrap1d(y, 1, -(y&1))*Lx + wrap1d(z, 2, -(z&1))*Lx*Ly;
  odd = wrap1d(x, 0, -((x&1)^1)) + wrap1d(y, 1, -((y&1)^1))*Lx + wrap1d(z, 2, -((z&1)^1))*Lx*Ly;
}

//Calculate the energy change for a rotation move
double Simulation::rotChange(int pos, int q) const {
  //Look at how many orientations match in old and new config
  int o = numOfNNOr(pos, array[pos].getOr());
  int n = numOfNNOr(pos, q);

  //de = -A*(n - o)
  return (double)(o - n)*A;
}

//Calculate the energy change for a swap move
double Simulation::swapChange(int pos1, int pos2) const {
  //Get id and orientation info
  int i1 = array[pos1].getId(), i2 = array[pos2].getId();
  int o1 = array[pos1].getOr(), o2 = array[pos2].getOr();

  //Calculate difference in id matches
  int did = (numOfNNId(pos1, i2) + numOfNNId(pos2, i1)) - (numOfNNId(pos1, i1) + numOfNNId(pos2, i2));

  //Calculate difference in orientation matches
  int dor = (numOfNNOr(pos1, o2) + numOfNNOr(pos2, o1)) - (numOfNNOr(pos1, o1) + numOfNNOr(pos2, o2));

  //de = -K((n_id1 + n_id2) - (o_id1 + o_id2)) - A((n_or1 + n_or2) - (o_or1 + o_or2))
  double de = -K*(double)did - A*(double)dor;

  //If the two are neighbors and are diff id or orient, need to subtract an overcount
  //This is because the environment was not actually updated
  if (areNN(pos1, pos2)) {
    if (i1 != i2) {
      de += 2*K;
    }
    if (o1 != o2) {
      de += 2*A;
    }
  }

  return de;
}

double Simulation::performMove(const Move* const m) const {
  int l0, l1, l2, l3;
  double e = 0.0;
  /* m->printMove(); */
  if (m->getType() == 0) {
    posLocks(m->getPos(), l0, l1);
    array[l0].lock();
    array[l1].lock();
    //Get energy change associated with rotation
    double de = rotChange(m->getPos(), m->getPar())/kT;

    //Check acceptance
    if (m->getProb() < exp(-de)) {
      //Update orientation and history
      array[m->getPos()].setOr(m->getPar());

      //Return change in energy
      e = de;
    }

    array[l1].unlock();
    array[l0].unlock();
  } else if (m->getType() == 1) {
    posLocks(m->getPos(), l0, l1);
    posLocks(m->getPar(), l2, l3);

    // Lock evens first, small then large
    if (l0 < l2) {
      array[l0].lock();
      array[l2].lock();
    } else if (l2 < l0) {
      array[l2].lock();
      array[l0].lock();
    }
    // Lock odds, small then large
    if (l1 < l3) {
      array[l1].lock();
      array[l3].lock();
    } else if (l3 < l1) {
      array[l3].lock();
      array[l1].lock();
    }

    //Calculate energy change associated with swap
    double de = swapChange(m->getPos(), m->getPar())/kT;

    //Check acceptance
    if (m->getProb() < exp(-de)) {
      //Swap cells
      swapIdOr(array[m->getPos()], array[m->getPar()]);

      //Return change in energy
      e = de;
    }

    // Unlock odds first, large then small
    if (l1 < l3) {
      array[l3].unlock();
      array[l1].unlock();
    } else if (l3 < l1) {
      array[l1].unlock();
      array[l3].unlock();
    }
    // Unlock evens, large then small
    if (l0 < l2) {
      array[l2].unlock();
      array[l0].unlock();
    } else if (l2 < l0) {
      array[l0].unlock();
      array[l2].unlock();
    }
  }
  return e;
}

//Calculate theta (M=26)
//From V. Argawal and B. Peters, J. Chem. Phys. 140, 084111
double* Simulation::calctheta() const {
  //Initialize
  double* theta = new double [NMAX];

  //Run through and calculate for every lattice position
  //TODO: parallel_for
  for (int i = 0; i < NMAX; i++) {
    theta[i] = 0.5;

    //Split into 3D coords
    int x = (i%(Lx*Ly))%Lx;
    int y = (i%(Lx*Ly))/Lx;
    int z = i/(Lx*Ly);

    //Run through neighbors (3x3x3 cube)
    int id = array[i].getId();
    for (int a = -1; a < 2; a++) {
      for (int b = -1; b < 2; b++) {
        for (int c = -1; c < 2; c++) {
          //Convert back to 1D coord
          int index = wrap1d(x, 0, a) + wrap1d(y, 1, b)*Lx + wrap1d(z, 2, c)*Lx*Ly;

          //Increase theta if id match
          if (id == array[index].getId()) {
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
double* Simulation::calcTheta() const {
  //Initialize
  double* Theta = new double [NMAX];

  //Get theta
  double* theta = calctheta();

  //TODO: parallel_for
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

//1D periodic boundary condition
int Simulation::wrap1d(int coord, int dir, int step) const {
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

//Calculate the energy between two lattice positions
double Simulation::pairEnergy(int pos1, int pos2) const {
  double e = 0.0;

  if (array[pos1].getId() == array[pos2].getId()) e -= K;
  if (array[pos1].getOr() == array[pos2].getOr()) e -= A;

  return e;
}

//Check to see how many neighbors of array[pos] have id i
int Simulation::numOfNNId(int pos, int i) const {
  int count = 0;

  //Run through neighbors and count matching ids
  for (int j = 0; j < 3; j++) {
    for (int k = -1; k <= 1; k += 2) {
      if (array[step3d(pos, j, k)].getId() == i) count++;
    }
  }

  return count;
}

//Check to see how many neighbors of array[pos] have orientation o
int Simulation::numOfNNOr(int pos, int o) const {
  int count = 0;

  //Run through neighbors and count matching ids
  for (int i = 0; i < 3; i++) {
    for (int j = -1; j <= 1; j += 2) {
      if (array[step3d(pos, i, j)].getOr() == o) count++;
    }
  }

  return count;
}

//Check to see if two indices are neighbors in 3D
int Simulation::areNN(int pos1, int pos2) const {
  for (int i = 0; i < 3; i++) {
    for (int j = -1; j <= 1; j+= 2) {
      if (step3d(pos1, i, j) == pos2) return 1;
    }
  }

  return 0;
}

//Swaps the identity and orientation of two cells
//This is cheaper than changing pointers and rearranging neighbor connections
void Simulation::swapIdOr(Cell& c1, Cell& c2) const {
  //Save 1's info
  int tempId = c1.getId(), tempOr = c1.getOr();

  //Pull 2's info into 1
  c1.setId(c2.getId()); c1.setOr(c2.getOr());

  //Put 1's info into 2
  c2.setId(tempId); c2.setOr(tempOr);

  return;
}

//Step of 1D index in 3D coords, with periodic boundary conditions
int Simulation::step3d(int index, int dir, int step) const {
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
Simulation::Simulation(int x, int y, int z, double T, double compA, double c) {
  Lx = x; Ly = y; Lz = z;
  NMAX = Lx*Ly*Lz;
  kT = T;
  cutoff = c;

  cout << "Hamiltonian: K=" << K << " A=" << A << endl;
  cout << "Move probabilities: Rotation=" << ROTATION << " Particle Swap=" << PARTSWAP << endl;
  cout << "Lattice dimensions of " << Lx << "x" << Ly << "x" << Lz << endl;
  cout << "Temperature (kT) of " << kT << endl;

  //Allocate array of Cells
  array = new Cell[NMAX];

  //Initialize cells
  int threshold = compA*NMAX;
  cout << "Simulation seeded with " << threshold << "/" << NMAX << " particles of species 1" << endl;
  cout << "Using Theta cutoff value of " << cutoff << endl;

  //Set species 1
  //TODO: parallel_for
  for (int i = 0; i < threshold; i++) {
    array[i].setId(1);
    array[i].setOr(1);
  }
  //Set species 2
  //TODO: parallel_for
  for (int i = threshold; i < NMAX; i++) {
    array[i].setId(2);
    array[i].setOr(1);
  }

  //Initialize random number generation
  auto seed = chrono::high_resolution_clock::now().time_since_epoch().count();
  cout << "Using seed " << seed << " for srand()" << endl;
  srand(seed);

  //Initialize energy of the system
  //TODO: Parallelize with atomic<double>?
  energy = 0.0;
  for (int i = 0; i < NMAX; i++) {
    //Calculate pairwise energy with forward pairs in each direction
    for (int j = 0; j < 3; j++) {
      energy += pairEnergy(i, step3d(i, j, 1));
    }
  }

  cout << "Initial energy: " << energy/kT << endl;
  cout << endl;
}

//Destructor
Simulation::~Simulation() {
  //Memory Management
  delete[] array;
}

//Function to evolve simulation by one sweep
void Simulation::doSweep() {
  /* uniform_int_distribution<> type(0, 1); */
  /* uniform_int_distribution<> pos(0, NMAX-1); */
  /* uniform_int_distribution<> orient(1, 6); */
  /* uniform_real_distribution<double> prob(0, 1); */
  const Move* moves[NMAX];
  for (int i = 0; i < NMAX; i++) {
    int ty = ((double) rand()/(double)RAND_MAX < ROTATION) ? 0 : 1;
    int p = rand()%NMAX;
    int q = (ty ? rand()%NMAX : rand()%6 + 1);
    double r = (double)rand()/(double)RAND_MAX;
    /* int ty = type(gen); */
    /* int p = pos(gen); */
    /* int q = (ty ? pos(gen) : orient(gen)); */
    /* double r = prob(gen); */
    moves[i] = new Move(ty, p, q, r);
    /* moves[i]->printMove(); */
  }
  energy += tbb::parallel_reduce(
      tbb::blocked_range<int>(0, NMAX),
      0.d,
      [&] (const tbb::blocked_range<int>& r, double init)->double {
      for (int i=r.begin(); i!=r.end(); ++i) {
        init += performMove(moves[i]);
        delete moves[i];
      }
      return init;
      }, [](double x, double y)->double {
      return x + y;
      });
  /* double de = 0.0; */
  /* #pragma omp parallel for reduction(+:de) */
  /* for (int i = 0; i < NMAX; i++) { */
  /*   de += performMove(moves[i]); */
  /*   delete moves[i]; */
  /* } */
  /* energy += de; */
}

//Function to return energy
double Simulation::getEnergy() const {
  return energy;
}

double* Simulation::calcThetaHistogram() const {
  //Initialize histogram
  double* h = new double [100];
  for (int i = 0; i < 100; i++) {
    h[i] = 0.0;
  }

  double* Theta = calcTheta();

  //TODO: Parallelize with atomic<double>?
  for (int i = 0; i < NMAX; i++) {
    //Get index in histogram
    int index = floorf(Theta[i]*100);
    if (index < 0) index = 0;
    if (index > 99) index = 99;

    h[index] += 1;
  }

  //Normalize
  //TODO: parallel_for
  for (int i = 0; i < 100; i++) {
    h[i] /= (double)NMAX;
  }

  //Memory Management
  delete[] Theta;

  return h;
}

double* Simulation::calcX1() const {
  int n [2] = {0, 0}; //Number of particles in each phase
  int n1 [2] = {0, 0}; // Number of species i in each phase

  double* Theta = calcTheta();

  //TODO: Parallelize with atomic<int>?
  for (int i = 0; i < NMAX; i++) {
    //Decide if in 1-rich or 2-rich phase
    if (Theta[i] >= cutoff) {
      n[0]++;
      if (array[i].getId() == 1) {
        n1[0]++;
      }
    } else {
      n[1]++;
      if (array[i].getId() == 1) {
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
