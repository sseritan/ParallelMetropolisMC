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
//Math include
#include <cmath>
//Vector include
#include <vector>
//Random include
#include <cstdlib>
//TBB
#include <tbb/parallel_for.h>

//Header include
#include "MetropolisMC.hpp"

using namespace std;

// Generate moves into M
void gen_moves(SimArray<int>& Y, vector<Move>* M, int& moves_left) {
  /* cout << "gen_moves.." << endl; */
  // (Re)initilize Y
  for(int i=0; i<Lx; i++) {
    for(int j=0; j<Ly; j++) {
      for(int k=0; k<Lz; k++) {
        // Initialize to 0
        Y[i][j][k] = 0;
      }
    }
  }

  // pc is used to fill moves in a round robin manner
  int pc = 0;
  while (true) {
    if (moves_left == 0) break;
    /* cout << "pc: " << pc << endl; */
    Move* move = new Move;

    move->cell0[0] = rand()%Lx;
    move->cell0[1] = rand()%Ly;
    move->cell0[2] = rand()%Lz;

    //Pick random number between 0 and 1 to determine move
    double m = (double)rand()/(double)RAND_MAX;

    //Choose move depending on m
    if (m < ROTATION) {
      //Rotation move
      move->type = 0;
    } else if (m > ROTATION && m < (ROTATION+PARTSWAP)) {
      move->type = 1;
      //Pick random particle to swap with
      move->cell1[0] = rand()%Lx;
      move->cell1[1] = rand()%Ly;
      move->cell1[2] = rand()%Lz;
    }

    int r = check_conflicts(Y, move, pc);
    /* cout << "check_conflicts returned " << r << endl; */
    if (r == 0) {
      // no conflict, assign round robin
      M[pc].push_back(*move);
      pc = (pc + 1) % np;
      moves_left--;
    } else if (r > 0) {
      // one conflict, assign to that processor
      M[r].push_back(*move);
      moves_left--;
    } else {
      // more than one conflict, we need to stop generate
      break;
    }
  }
}

void print_move(Move& move) {
  if (move.type) { // swap
    cout << "swap((" << move.cell0[0] << ", " << move.cell0[1] << ", " <<
      move.cell0[2] << "), (" << move.cell1[0] << ", " << move.cell1[1] << ", "
      << move.cell1[2] << "))" << endl;
  } else { // rotation
    cout << "rotate(" << move.cell0[0] << ", " << move.cell0[1] << ", " <<
      move.cell0[2] << ")" << endl;
  }
}

// Check wheter a move conflicts, return what processor a move should be assigned to
// return -1 if there are two or more conflicts
int check_conflicts(SimArray<int>& Y, Move* move, int pc) {
  /* cout << "check conflicts for move:" << endl; */
  /* print_move(*move); */
  int cell = 0;
  int conflict_processor = 0;
  if (move->type) { // swap
    for (int i = 0; i < 3; i++) {
      for (int j = -1; j <= 1; j += 2) {
        cell = Y[mod(move->cell0[0] + (i==2)*j, Lx)][mod(move->cell0[1] +
            (i==1)*j, Ly)][mod(move->cell0[2] + (i==0)*j, Lz)];
        /* cout << "cell coordinates: " << mod(move->cell0[0] + (i==2)*j, Lx) << */
        /*   ", " << mod(move->cell0[1] + (i==1)*j, Ly) << ", " << */
        /*   mod(move->cell0[2] + (i==0)*j, Lz) << endl; */
        /* cout << "cell content: " << cell << endl; */
        if (cell) {
          if (conflict_processor and (conflict_processor != cell)) {
            return -1; // more than one conflicts!
          } else {
            conflict_processor = cell;
          }
        }
      }
    }
    for (int i = 0; i < 3; i++) {
      for (int j = -1; j <= 1; j += 2) {
        cell = Y[mod(move->cell1[0] + (i==2)*j, Lx)][mod(move->cell1[1] +
            (i==1)*j, Ly)][mod(move->cell1[2] + (i==0)*j, Lz)];
        /* cout << "cell coordinates: " << mod(move->cell1[0] + (i==2)*j, Lx) << */
        /*   ", " << mod(move->cell1[1] + (i==1)*j, Ly) << ", " << */
        /*   mod(move->cell1[2] + (i==0)*j, Lz) << endl; */
        /* cout << "cell content: " << cell << endl; */
        if (cell) {
          if (conflict_processor and (conflict_processor != cell)) {
            return -1; // more than one conflicts!
          } else {
            conflict_processor = cell;
          }
        }
      }
    }
    if (conflict_processor > 0) pc = conflict_processor;
    for (int i = 0; i < 3; i++) {
      for (int j = -1; j <= 1; j += 2) {
        Y[mod(move->cell0[0] + (i==2)*j, Lx)][mod(move->cell0[1] +
            (i==1)*j, Ly)][mod(move->cell0[2] + (i==0)*j, Lz)] = pc;
        /* cout << "set cell to " << pc << endl; */
      }
    }
    for (int i = 0; i < 3; i++) {
      for (int j = -1; j <= 1; j += 2) {
        Y[mod(move->cell1[0] + (i==2)*j, Lx)][mod(move->cell1[1] +
            (i==1)*j, Ly)][mod(move->cell1[2] + (i==0)*j, Lz)] = pc;
        /* cout << "set cell to " << pc << endl; */
      }
    }
  } else { // rotation
    for (int i = 0; i < 3; i++) {
      for (int j = -1; j <= 1; j += 2) {
        cell = Y[mod(move->cell0[0] + (i==2)*j, Lx)][mod(move->cell0[1] +
            (i==1)*j, Ly)][mod(move->cell0[2] + (i==0)*j, Lz)];
        /* cout << "cell coordinates: " << mod(move->cell0[0] + (i==2)*j, Lx) << */
        /*   ", " << mod(move->cell0[1] + (i==1)*j, Ly) << ", " << */
        /*   mod(move->cell0[2] + (i==0)*j, Lz) << endl; */
        /* cout << "cell content: " << cell << endl; */
        if (cell) {
          if (conflict_processor and (conflict_processor != cell)) {
            return -1; // more than one conflicts!
          } else {
            conflict_processor = cell;
          }
        }
      }
    }
    if (conflict_processor > 0) pc = conflict_processor;
    for (int i = 0; i < 3; i++) {
      for (int j = -1; j <= 1; j += 2) {
        Y[mod(move->cell0[0] + (i==2)*j, Lx)][mod(move->cell0[1] +
            (i==1)*j, Ly)][mod(move->cell0[2] + (i==0)*j, Lz)] = pc;
        /* cout << "set cell to " << pc << endl; */
      }
    }
  }
  return conflict_processor;
}

//Run full sweep (Lx*Ly*Lz MC moves)
void sweep(vector<Move>* M, SimArray<int>& X, SimArray<int>& S, double& e, const double kT) {
  /* cout << "sweep" << endl; */
  int m = Lx*Ly*Lz;
  int count = 0;

  // Y is used to keep track of move assignment
  SimArray<int>* Y_ptr = new SimArray<int>;

  while (m > 0) { // until we do enough moves
    count++;

    gen_moves(*Y_ptr, M, m);

    /* cout << "moves left for this sweep: " << m << endl; */

    tbb::parallel_for (0, np, [&] (int i) {
      do_moves(M[i], X, S, e, kT);
    });
  }
  /* cout << "average batch size: " << Lx*Ly*Lz / count << endl; */
  delete Y_ptr;
}

void do_moves(vector<Move>& Mi, SimArray<int>& X, SimArray<int>& S, double& e, const double kT) {
  /* cout << "making moves for i = " << i << endl; */
  while (Mi.size() > 0) {
    Move move = Mi.back();
    Mi.pop_back();
    /* print_move(move); */
    if (move.type) { // swap
      //Get change in energy associated with switching particles
      int i = move.cell0[0];
      int j = move.cell0[1];
      int k = move.cell0[2];
      int ii = move.cell1[0];
      int jj = move.cell1[1];
      int kk = move.cell1[2];
      double de = particle_swap_energy_change(X, S, i, j, k, ii, jj, kk)/kT;

      //Check whether to accept move or not
      if ((double)rand()/(double)RAND_MAX < exp(-de)) {
        //Accept move. Update orientation and energy
        int x = X[i][j][k]; int s = S[i][j][k];
        X[i][j][k] = X[ii][jj][kk]; S[i][j][k] = S[ii][jj][kk];
        X[ii][jj][kk] = x; S[ii][jj][kk] = s;
        e += de;
      }
    } else { // rotation
      int i = move.cell0[0];
      int j = move.cell0[1];
      int k = move.cell0[2];
      //Pick new orientation
      int q = rand()%6 + 1;

      //Calculate energy change for chosen rotation
      double de = rotation_energy_change(S, i, j, k, q)/kT;

      //Check whether to accept or not
      if ((double)rand()/(double)RAND_MAX < exp(-de)) {
        //Accept move. Update orientation and energy
        S[i][j][k] = q;
        e += de;
      }
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

//Calculate the energy between two particles (e = -K dm(i)m(j) -A ds(i)s(j)
double pairwise_energy(int m1, int m2, int s1, int s2) {
  double e = 0.0;

  //Identity bonus
  if (m1 == m2) e += -K;

  //Orientation bonus
  if (s1 == s2) e += -A;

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

//Calculate the energy in a rotation move
double rotation_energy_change(const SimArray<int>& S, int i, int j, int k, int q2) {
  //Store original orientation
  int q1 = S[i][j][k];

  int qOrig = 0, qNew = 0;
  //Only calculate if it actually makes a difference
  if (q1 != q2) {
    //Store neighboring orientations for ease
    vector<int> o;
    o.push_back(S[mod(i-1,Lx)][j][k]);
    o.push_back(S[mod(i+1,Lx)][j][k]);
    o.push_back(S[i][mod(j-1,Ly)][k]);
    o.push_back(S[i][mod(j+1,Ly)][k]);
    o.push_back(S[i][j][mod(k-1,Lz)]);
    o.push_back(S[i][j][mod(k+1,Lz)]);

    //Count number of orientation matches in both cases
    for (int q : o) {
      if (q == q1) {
        qOrig++;
      }
      if (q == q2) {
        qNew++;
      }
    }
  }

  //Energy change = -A(qNew - qOrig)
  return A*(qOrig - qNew);
}

//Calculate energy in a particle swap move
//Note that this is not the most efficient, but not bottleneck
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

//Calculate full arrays for phase parameters (TODO: Highly parallelizable)
SimArray<double> phase_parameter(const SimArray<int>& X) {
  //Initialize temp arrays
  SimArray<double> theta, Theta;

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
double local_phase_parameter(const SimArray<int>& X, int i, int j, int k) {
  //Initialize parameter
  double theta = 0.5;

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
double medium_phase_parameter(const SimArray<double>& theta, int i, int j, int k) {
  //Initialize parameter
  double Theta = 0.0;

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

//From the given respective CUTOFF value, calculate XA in the two phases (0 = A-rich, 1 = B-rich)
array<double, 2> phase_data(const SimArray<int>& X, const SimArray<double>& param, const double cutoff) {
  //Temporary number arrays (total and A in both phases
  array<int, 2> n = {0, 0};
  array<int, 2> nA = {0, 0};

  //Count and get numbers for both phases
  for (int i = 0; i < Lx; i++) {
    for (int j = 0; j < Ly; j++) {
      for (int k = 0; k < Lz; k++) {
        if (param[i][j][k] >= cutoff) {
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
  array<double, 2> XA;
  for (int i = 0; i < 2; i++) {
    XA[i] = (double)nA[i]/(double)n[i];
  }

  return XA;
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
