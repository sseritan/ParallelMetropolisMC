//
//
// Originally written by Stefan Seritan on 8/28/14
// Modified by Stefan Seritan and Wei Dai for CS 140 Final Project Winter 2015
//
// All subroutines required for MetropolisMC.cpp
//

//I/O include
#include <iostream>
#include <fstream>
//Array include
#include <array>
//Stdlib include for srand
#include <cstdlib>
//Math include
#include <cmath>

//Header include
#include "MetropolisMC.hpp"

using namespace std;

/* FUNCTION PROTOTYPES */
void sweep(SimArray<int>& X, SimArray<int>& S, double& e);

double energy(const SimArray<int>& X, const SimArray<int>& S);
double pairwise_energy(int m1, int m2, int s1, int s2);
double point_energy(const SimArray<int>& X, const SimArray<int>& S, int i, int j, int k);
double rotation_energy_change(const SimArray<int>& X, SimArray<int>& S, int i, int j, int k, int q);
double particle_swap_energy_change(SimArray<int>& X, SimArray<int>& S, int i, int j, int k, int ii, int jj, int kk);

SimArray<float> phase_parameter(const SimArray<int>& X);
float local_phase_parameter(const SimArray<int>& X, int i, int j, int k);
float medium_phase_parameter(const SimArray<float>& theta, int i, int j, int k);
SimArray<float> orientation_parameter(const SimArray<int>& S);
float local_orientation_parameter(const SimArray<int>& S, int i, int j, int k);
float medium_orientation_parameter(const SimArray<float>& phi, int i, int j, int k);

array<float, 100> histogram(const SimArray<float>& param);
Dim2Array histogram2d(const SimArray<float>& Theta, const SimArray<float>& Phi);
array<float, 2> phase_data(const SimArray<int>& X, const SimArray<float>& Theta);

int mc_acc(const double de);
int mod(int k, int n);
int random_int(int a, int b);
double rand01();

void print_VMD_snapshot(SimArray<int>& X, SimArray<int>& S, int t);

/* FUNCTION DEFINITIONS */
//Run full sweep (Lx*Ly*Lz MC moves)
void sweep(SimArray<int>& X, SimArray<int>& S, double& e) {
  for (int t=0; t < Lx*Ly*Lz; t++) {
    //Pick random lattice point
    int i = random_int(0, Lx-1);
    int j = random_int(0, Ly-1);
    int k = random_int(0, Lz-1);

    if (DEBUGGING) {
      cout << "i:" << i << " j:" << j << " k:" << k << " ";
    }

    //Pick random number between 0 and 1 to determine move
    double m = rand01();

    //Choose move depending on m
    if (m < ROTATION) {
      //Rotation move
      //Pick new orientation
      int q = random_int(1, Q);

      //Calculate energy change for chosen rotation
      double de = rotation_energy_change(X, S, i, j, k, q)/kT;

      //Get acceptance
      int acc = mc_acc(de);

      //Print out move information if DEBUGGING flag set
      if (DEBUGGING) {
        cout << "Rot: qOrg=" << S[i][j][k] << " q=" << q << " de=" << de << (acc ? " Accepted" : " Rejected");
      }

      //Check whether to accept or not
      if (acc) {
        //Accept move. Update orientation and energy
        S[i][j][k] = q;
        e += de;
      }
    } else if (m > ROTATION && m < (ROTATION+PARTSWAP)) {
      //Pick random particle to swap with
      int ii = random_int(0, Lx-1);
      int jj = random_int(0, Ly-1);
      int kk = random_int(0, Lz-1);

      //Get change in energy associated with switching particles
      double de = particle_swap_energy_change(X, S, i, j, k, ii, jj, kk)/kT;

      //Get acceptance
      int acc = mc_acc(de);

      //Print out move information if DEBUGGING flag set
      if (DEBUGGING) {
        cout << " ii:" << ii << " jj:" << jj << " kk:" << kk;
        cout << " PartSwap: m1=" << X[i][j][k] << " m2=" << X[ii][jj][kk] << " s1=" << S[i][j][k] <<
          " s2=" << S[ii][jj][kk] << " de=" << de << (acc ? " Accepted" : " Rejected");
      }

      //Check whether to accept move or not
      if (acc) {
        //Accept move. Update orientation and energy
        int x = X[i][j][k]; int s = S[i][j][k];
        X[i][j][k] = X[ii][jj][kk]; S[i][j][k] = S[ii][jj][kk];
        X[ii][jj][kk] = x; S[ii][jj][kk] = s;
        e += de;
      }
    }

    if (DEBUGGING) {
      cout << endl;
    }
  }
}

//Calculate energy of full box
double energy(const SimArray<int>& X, const SimArray<int>& S) {
  double e = 0.0;

  //Iterate through box (with periodic boundaries)
  for(int i=0; i<Lx; i++) {
    int ii = mod(i+1, Lx);
    for(int j=0; j<Ly; j++) {
      int jj = mod(j+1, Ly);
      for(int k=0; k<Lz; k++) {
        int kk = mod(k+1, Lz);

        //Calculate pairwise energy (forward in each direction)
        e += pairwise_energy(X[i][j][k], X[ii][j][k], S[i][j][k], S[ii][j][k]);
        e += pairwise_energy(X[i][j][k], X[i][jj][k], S[i][j][k], S[i][jj][k]);
        e += pairwise_energy(X[i][j][k], X[i][j][kk], S[i][j][k], S[i][j][kk]);
      }
    }
  }

  return e;
}

//Calculate the energy between two particles
double pairwise_energy(int m1, int m2, int s1, int s2) {
  double e = 0.0;

  //Hamiltonian is -K (if id eq) - A (if orient eq)

  //Identity bonus
  if (m1 == m2) {
    e += -K;
  }

  //Orientation bonus
  if (s1 == s2) {
    e += -A;
  }

  return e;
}

//Calculate the energy with nearest neighbors at point (i,j,k)
double point_energy(const SimArray<int>& X, const SimArray<int>& S, int i, int j, int k) {
  //Calculate each pairwise energy
  double e = pairwise_energy(X[i][j][k], X[mod(i+1, Lx)][j][k], S[i][j][k], S[mod(i+1, Lx)][j][k]);
  e += pairwise_energy(X[i][j][k], X[mod(i-1, Lx)][j][k], S[i][j][k], S[mod(i-1, Lx)][j][k]);
  e += pairwise_energy(X[i][j][k], X[i][mod(j+1, Ly)][k], S[i][j][k], S[i][mod(j+1, Ly)][k]);
  e += pairwise_energy(X[i][j][k], X[i][mod(j-1, Ly)][k], S[i][j][k], S[i][mod(j-1, Ly)][k]);
  e += pairwise_energy(X[i][j][k], X[i][j][mod(k+1, Lz)], S[i][j][k], S[i][j][mod(k+1, Lz)]);
  e += pairwise_energy(X[i][j][k], X[i][j][mod(k-1, Lz)], S[i][j][k], S[i][j][mod(k-1, Lz)]);

  return e;
}

//Calculate the energy change for assuming the orientation q at position (i,j,k)
double rotation_energy_change(const SimArray<int>& X, SimArray<int>& S, int i, int j, int k, int q) {
  //Calculate energy with original orientation (with 6 nearest neighbors)
  double e1 = point_energy(X, S, i, j, k);

  //Save original orientation, and update with new orientation
  int qOrg = S[i][j][k];
  S[i][j][k] = q;

  //Calculate energy with new orientation
  double e2 = point_energy(X, S, i, j, k);

  //Restore old orientation
  S[i][j][k] = qOrg;

  return (e2-e1);
}

//Calculate the energy change for switching particles at positions (i,j,k) and (ii,jj,kk)
double particle_swap_energy_change(SimArray<int>& X, SimArray<int>& S, int i, int j, int k, int ii, int jj, int kk) {
  //Calculate energy with original orientation
  double e1 = point_energy(X, S, i, j, k) + point_energy(X, S, ii, jj, kk);

  //Save original information
  int m1 = X[i][j][k], m2 = X[ii][jj][kk];
  int s1 = S[i][j][k], s2 = S[ii][jj][kk];

  //Update identity and orientation
  X[i][j][k] = m2; X[ii][jj][kk] = m1;
  S[i][j][k] = s2; S[ii][jj][kk] = s1;

  //Calculate energy with new orientation
  double e2 = point_energy(X, S, i, j, k) + point_energy(X, S, ii, jj, kk);

  //Restore old identity and orientation
  X[i][j][k] = m1; X[ii][jj][kk] = m2;
  S[i][j][k] = s1; S[ii][jj][kk] = s2;

  return (e2-e1);
}

//Calculate full arrays for phase parameters (TODO: Highly parallelizable) TODO: Fix return by value?
SimArray<float> phase_parameter(const SimArray<int>& X) {
  //Initialize temp arrays
  SimArray<float> theta, Theta;

  //Calculate theta from x
  for (int i = 0; i < Lx; i++) {
    for (int j = 0; j < Ly; j++) {
      for (int k = 0; k < Lz; k++) {
        theta[i][j][k] = local_phase_parameter(X, i, j, k);
      }
    }
  }

  //Calculate Theta from theta
  for (int i = 0; i < Lx; i++) {
    for (int j = 0; j < Ly; j++) {
      for (int k = 0; k < Lz; k++) {
        Theta[i][j][k] = medium_phase_parameter(theta, i, j, k);
      }
    }
  }

  return Theta;
}

//Calculate theta, the local order parameter for determining phase
//From V. Argawal and B. Peters, J. Chem. Phys. 140, 084111
float local_phase_parameter(const SimArray<int>& X, int i, int j, int k) {
  //Initialize parameter
  float theta = 0.5;

  //Sum over nearest, next-nearest, and next-next-neighbors (3x3x3 cube)
  int id = X[i][j][k];
  for (int ip = -1; ip <= 1; ip++) {
    for (int jp = -1; jp <= 1; jp++) {
      for (int kp = -1; kp <= 1; kp++) {
        //Add 1/(2*26) for same species, -1/(2*26) for opposite species
        if (id == X[mod(i+ip, Lx)][mod(j+jp, Ly)][mod(k+kp, Lz)]) {
          if (id == 1) {
            theta += 1.0/52.0;
          } else {
            theta += -1.0/52.0;
          }
        }
      }
    }
  }

  //Remove overcounting of itself
  if (id == 1) {
    theta += -1.0/52.0;
  } else {
    theta += 1.0/52.0;
  }

  //Sometimes sign is weird, gives -0.00
  if (theta < 0.0) {
    theta = 0.0;
  }

  return theta;
}

//Calculate Theta, the medium order parameter for determining phase (average of theta)
//From V. Agarwal and B. Peters, J. Chem. Phys. 140, 084111
float medium_phase_parameter(const SimArray<float>& theta, int i, int j, int k) {
  //Initialize parameter
  float Theta = 0.0;

  //Run through array and take average over 3x3x3 cube
  for (int ip = -1; ip <= 1; ip++) {
    for (int jp = -1; jp <= 1; jp++) {
      for (int kp = -1; kp <= 1; kp++) {
        Theta += theta[mod(i+ip, Lx)][mod(j+jp, Ly)][mod(k+kp, Lz)];
      }
    }
  }

  //Average out Theta
  return Theta/27.0;
}

//Calculate full arrays of the orientation parameters //TODO:Fix return by value?
SimArray<float> orientation_parameter(const SimArray<int>& S) {
  //Initialize temp arrays
  SimArray<float> phi, Phi;

  //Calculate phi from s
  for (int i = 0; i < Lx; i++) {
    for (int j = 0; j < Ly; j++) {
      for (int k = 0; k < Lz; k++) {
        phi[i][j][k] = local_orientation_parameter(S, i, j, k);
      }
    }
  }

  //Calculate Phi from phi
  for (int i = 0; i < Lx; i++) {
    for (int j = 0; j < Ly; j++) {
      for (int k = 0; k < Lz; k++) {
        Phi[i][j][k] = medium_orientation_parameter(phi, i, j, k);
      }
    }
  }

  return Phi;
}

//Calculate the local orientation parameter (based on local phase parameter)
float local_orientation_parameter(const SimArray<int>& S, int i, int j, int k) {
  //Initialize variable
  float phi = 0.0;

  //Save orientation of (i,j,k)
  int q = S[i][j][k];
  for (int a = -1; a <= 1; a++) {
    for (int b = -1; b <= 1; b++) {
      for (int c = -1; c <= 1; c++) {
        if (q == S[mod(i+a, Lx)][mod(j+b, Ly)][mod(k+c, Lz)]) {
          phi += 1.0;
        }
      }
    }
  }

  //Remove overcounting of center
  phi += -1.0;

  //Average out phi
  return phi/26.0;
}

//Calculate the medium orientation parameter (based on medium phase parameter)
float medium_orientation_parameter(const SimArray<float>& phi, int i, int j, int k) {
  //Initialize variable
  float Phi = 0.0;

  //Run through phi array
  for (int a = -1; a <= 1; a++) {
    for (int b = -1; b <= 1; b++) {
      for (int c = -1; c <= 1; c++) {
        Phi += phi[i][j][k];
      }
    }
  }

  //Average out Phi
  return Phi/27.0;
}

//Make a histogram of an order parameter
array<float, 100> histogram(const SimArray<float>& param) {
  //Initialize array to all zeroes
  array<float, 100> h = {0.0};

  //Run through parameter array
  for (int i = 0; i < Lx; i++) {
    for (int j = 0; j < Ly; j++) {
      for (int k = 0; k < Lz; k++) {
        //Truncate parameter at 2nd decimal place to get bin
        int index = floorf(param[i][j][k]*100);

        //1.00 will return 100, but 99 is max bin
        if (index > 99) {
          index = 99;
        }

        //Add 1 to corresponding bin
        h[index] += 1;
      }
    }
  }

  //Normalize histogram
  for (int i = 0; i < 100; i++) {
    h[i] = h[i]/(float)(Lx*Ly*Lz);
  }

  return h;
}

Dim2Array histogram2d(const SimArray<float>& Theta, const SimArray<float>& Phi) {
  //Initialize 2d histogram
  Dim2Array h = {{0.0}};

  for (int i = 0; i < Lx; i++) {
    for (int j = 0; j < Ly; j++) {
      for (int k = 0; k < Lz; k++) {
        int ThetaIndex = floorf(Theta[i][j][k]*100);
        int PhiIndex = floorf(Phi[i][j][k]*100);

        //Make sure not out of bounds
        if (ThetaIndex > 99) {
          ThetaIndex = 99;
        }
        if (PhiIndex > 99) {
          PhiIndex = 99;
        }

        //Increment histogram. Theta is x, so inner level
        h[PhiIndex][ThetaIndex] += 1.0;
      }
    }
  }

  //Normalize histogram
  for (auto a : h) {
    for (float f : a) {
      f = f/(float)(Lx*Ly*Lz);
    }
  }

  return h;
}

//From the given respective CUTOFF value, calculate XA in the two phases (0 = A-rich, 1 = B-rich)
array<float, 2> phase_data(const SimArray<int>& X, const SimArray<float>& param) {
  //Temporary number arrays (total and A in both phases
  array<int, 2> n = {0, 0}; array<int, 2> nA = {0, 0};

  //Count and get numbers for both phases
  //Note the ternary operator lets me do both cases (by Theta or by Phi)
  //This is because large Theta and large Phi both correspond to the A-rich phase (given my Hamiltonian)
  for (int i = 0; i < Lx; i++) {
    for (int j = 0; j < Ly; j++) {
      for (int k = 0; k < Lz; k++) {
        if (param[i][j][k] >= ((RUNTYPE == 3) ? THETA_CUTOFF : PHI_CUTOFF)) {
          n[0]++;
          if (X[i][j][k] == 1) {
            nA[0]++;
          }
        } else {
          n[1]++;
          if (X[i][j][k] == 1) {
            nA[1]++;
          }
        }
      }
    }
  }

  //Calculate mole fraction in both phases
  array<float, 2> XA;
  for (int i = 0; i < 2; i++) {
    XA[i] = (float)nA[i]/(float)n[i];
  }

  return XA;
}

//Metropolis Monte Carlo acceptance criteria. Returns 1 if accepted, 0 if rejected.
int mc_acc(const double de) {
  //Get random cutoff number
  double r = rand01();

  //If probability of acceptance is smaller than machine precision, reject
  if (de*kT > 999.0) {
    return 0;
  }

  if (DEBUGGING) {
    cout << " r: " << r << " exp(-de): " << exp(-de);
  }

  //Accept if greater than r
  if (r < exp(-de)) {
    return 1;
  }

  return 0;
}

//Modified modulus function that guarantees periodic boundary conditions
int mod(int k, int n) {
  //Keep k between 0 and n-1
  while(k >= n) {
    k -= n;
  }
  while(k < 0) {
    k += n;
  }

  return k;
}

//Baron's random integer generator
int random_int(int a, int b) {
  //Return a random number in the given range [a,b]
  return (rand()%(b-a+1) + a);
}

//Baron's random double in interval [0,1] generator
double rand01() {
  int max = 1000000;
  int r = random_int(1, max);
  return ((double)r-0.5)/(double)max;
}

//Print out lattice to VMD_snapshot.txt
void print_VMD_snapshot(SimArray<int>& X, SimArray<int>& S, int t) {
  //Open file
  ofstream VMDFILE ("VMD_snapshot.txt", ios::app);

  //Timestep
  VMDFILE << "ITEM: TIMESTEP\n" << t << endl;

  //Number of atoms
  VMDFILE << "ITEM: NUMBER OF ATOMS\n" << Lx*Ly*Lz << endl;

  //Bounds
  VMDFILE << "ITEM: BOX BOUNDS pp pp pp" << endl;
  VMDFILE << "0 " << Lx << "\n0 " << Ly << "\n0 " << Lz << endl;

  //Lattice output
  VMDFILE << "ITEM: ATOMS id vx xs ys zs" << endl;
  for (int i = 0; i < Lx; i++) {
    for (int j = 0; j < Ly; j++) {
      for (int k = 0; k < Lz; k++) {
        VMDFILE << (Ly*Lz*i + Lz*j + k) << " " << (10*X[i][j][k]+S[i][j][k]) << " ";
        VMDFILE << i << " " << j << " " << k << endl;
      }
    }
  }

  //Close file
  VMDFILE.close();
}
