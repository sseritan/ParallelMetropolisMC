//
// Written by Stefan Seritan on 03/08/15 for CS 140
// Contains the SimTest class, used for testing Simulation class
//

#include <iostream>

#include "./CellMove.hpp"
#include "./Simulation.hpp"
#include "./SimTest.hpp"

using namespace std;

//Dummy constructor
SimTest::SimTest() : Simulation(4, 4, 4, 1.0, 1.0, 1.0) {
  //Hardcode some cells for testing (id and orient are the same for simplicity)
  //k = 0    k = 1    k = 2    k=3
  //1 1 1 1  1 1 2 1  1 2 2 2  1 1 2 1
  //1 1 1 2  1 1 2 1  1 2 2 2  1 1 2 1
  //2 1 2 1  1 2 1 1  2 1 2 1  1 2 1 1
  //2 1 1 2  1 1 1 1  1 2 2 1  2 1 1 1
  array[3 + 1*Lx + 0*Lx*Ly]->setId(2); array[3 + 1*Lx + 0*Lx*Ly]->setOr(2);
  array[0 + 2*Lx + 0*Lx*Ly]->setId(2); array[0 + 2*Lx + 0*Lx*Ly]->setOr(2);
  array[2 + 2*Lx + 0*Lx*Ly]->setId(2); array[2 + 2*Lx + 0*Lx*Ly]->setOr(2);
  array[0 + 3*Lx + 0*Lx*Ly]->setId(2); array[0 + 3*Lx + 0*Lx*Ly]->setOr(2);
  array[3 + 3*Lx + 0*Lx*Ly]->setId(2); array[3 + 3*Lx + 0*Lx*Ly]->setOr(2);

  array[2 + 0*Lx + 1*Lx*Ly]->setId(2); array[2 + 0*Lx + 1*Lx*Ly]->setOr(2);
  array[2 + 1*Lx + 1*Lx*Ly]->setId(2); array[2 + 1*Lx + 1*Lx*Ly]->setOr(2);
  array[1 + 2*Lx + 1*Lx*Ly]->setId(2); array[1 + 2*Lx + 1*Lx*Ly]->setOr(2);

  array[1 + 0*Lx + 2*Lx*Ly]->setId(2); array[1 + 0*Lx + 2*Lx*Ly]->setOr(2);
  array[2 + 0*Lx + 2*Lx*Ly]->setId(2); array[2 + 0*Lx + 2*Lx*Ly]->setOr(2);
  array[3 + 0*Lx + 2*Lx*Ly]->setId(2); array[3 + 0*Lx + 2*Lx*Ly]->setOr(2);
  array[1 + 1*Lx + 2*Lx*Ly]->setId(2); array[1 + 1*Lx + 2*Lx*Ly]->setOr(2);
  array[2 + 1*Lx + 2*Lx*Ly]->setId(2); array[2 + 1*Lx + 2*Lx*Ly]->setOr(2);
  array[3 + 1*Lx + 2*Lx*Ly]->setId(2); array[3 + 1*Lx + 2*Lx*Ly]->setOr(2);
  array[0 + 2*Lx + 2*Lx*Ly]->setId(2); array[0 + 2*Lx + 2*Lx*Ly]->setOr(2);
  array[2 + 2*Lx + 2*Lx*Ly]->setId(2); array[2 + 2*Lx + 2*Lx*Ly]->setOr(2);
  array[1 + 3*Lx + 2*Lx*Ly]->setId(2); array[1 + 3*Lx + 2*Lx*Ly]->setOr(2);
  array[2 + 3*Lx + 2*Lx*Ly]->setId(2); array[2 + 3*Lx + 2*Lx*Ly]->setOr(2);

  array[2 + 0*Lx + 3*Lx*Ly]->setId(2); array[2 + 0*Lx + 3*Lx*Ly]->setOr(2);
  array[2 + 1*Lx + 3*Lx*Ly]->setId(2); array[2 + 1*Lx + 3*Lx*Ly]->setOr(2);
  array[1 + 2*Lx + 3*Lx*Ly]->setId(2); array[1 + 2*Lx + 3*Lx*Ly]->setOr(2);
  array[0 + 3*Lx + 3*Lx*Ly]->setId(2); array[0 + 3*Lx + 3*Lx*Ly]->setOr(2);

  //Print out array just in case
  /*for (int k = 0; k < 4; k++) {
    for (int j = 0; j < 4; j++) {
      for (int i = 0; i < 4; i++) {
        cout << array[i+ Lx*j + Lx*Ly*k]->getId() << array[i+ Lx*j + Lx*Ly*k]->getOr() << " ";
      }
      cout << endl;
    }
    cout << endl;
  }*/

  //Initialize energies
  energy = 0.0;
  for (int i = 0; i < NMAX; i++) {
    //Calculate pairwise energy with forward pairs in each direction
    for (int j = 0; j < 3; j++) {
      energy += pairEnergy(i, step3d(i, j, 1));
    }
  }
}

int SimTest::compExp(double r, double e) {
  if (r != e) {
    cout << "Failed. Expected " << e << " but got " << r << "." << endl;
    return 1;
  }
  cout << "Passed!" << endl;
  return 0;
}

//Test wrap1d function
int SimTest::wrap() {
  int failed = 0;
  int f = 0;

  //Wrap in i direction
  cout << "wrap1d in 0 (i) direction: " << endl;
  for (int i = 0; i < Lx; i++) {
    f = compExp(wrap1d(i, 0, 1), ((i == (Lx-1)) ? 0 : i+1));
    if (f) cout << "Failed at " << i << " ";
    failed += f;
  }
  f = compExp(wrap1d(0, 0, -1), 3);
  if (f) cout << "Failed at 0 going backwards" << endl;
  failed += f;

  //Wrap in j direction
  cout << "wrap1d in 1 (j) direction: " << endl;
  for (int i = 0; i < Ly; i++) {
    f = compExp(wrap1d(i, 1, 1), ((i == Ly-1) ? 0 : i+1));
    if (f) cout << "Failed at " << i << " ";
    failed += f;
  }
  f = compExp(wrap1d(0, 1, -1), 3);
  if (f) cout << "Failed at 0 going backwards" << endl;
  failed += f;

  //Wrap in j direction
  cout << "wrap1d in 2 (k) direction: " << endl;
  for (int i = 0; i < Lz; i++) {
    f = compExp(wrap1d(i, 2, 1), ((i == Lz-1) ? 0 : i+1));
    if (f) cout << "Failed at " << i << " ";
    failed += f;
  }
  f = compExp(wrap1d(0, 2, -1), 3);
  if (f) cout << "Failed at 0 going backwards" << endl;
  failed += f;

  return failed;
}

//Test step3d function
int SimTest::step() {
  int failed = 0;

  cout << "i-1 (center): ";
  failed += compExp(step3d(5, 0, -1), 4);

  cout << "i+1 (center): ";
  failed += compExp(step3d(5, 0, 1), 6);

  cout << "j-1 (center): ";
  failed += compExp(step3d(5, 1, -1), 1);

  cout << "j+1 (center): ";
  failed += compExp(step3d(5, 1, 1), 9);

  cout << "k-1 (center): ";
  failed += compExp(step3d(5, 2, -1), 53);

  cout << "k+1 (center): ";
  failed += compExp(step3d(5, 2, 1), 21);

  cout << "i-1 (origin): ";
  failed += compExp(step3d(0, 0, -1), 3);

  cout << "i+1 (origin): ";
  failed += compExp(step3d(0, 0, 1), 1);

  cout << "j-1 (origin): ";
  failed += compExp(step3d(0, 1, -1), 12);

  cout << "j+1 (origin): ";
  failed += compExp(step3d(0, 1, 1), 4);

  cout << "k-1 (origin): ";
  failed += compExp(step3d(0, 2, -1), 48);

  cout << "k+1 (origin): ";
  failed += compExp(step3d(0, 2, 1), 16);

  return failed;
}

//Test numOfNNId function
int SimTest::nn_id_count() {
  int failed = 0;

  //(1, 1, 0) 1 in 1s
  cout << "numOfNNId (matching 1, central): ";
  failed += compExp(numOfNNId(5, 1), 6);

  //(1, 0, 0) 1 in 1s
  cout << "numOfNNId (matching 1, boundary): ";
  failed += compExp(numOfNNId(1, 1), 6);

  //(2, 1, 2) 2 in 2s
  cout << "numOfNNId (matching 2, central): ";
  failed += compExp(numOfNNId(38, 2), 6);

  //(2, 0, 2) 2 in 2s
  cout << "numOfNNId (matching 2, boundary): ";
  failed += compExp(numOfNNId(34, 2), 6);

  //(0, 3, 0) 2 in mixed
  cout << "numOfNNId (perfect mix, boundary): ";
  failed += compExp(numOfNNId(12, 2), 3);

  return failed;
}

//Test numOfNNOr function
int SimTest::nn_or_count() {
  int failed = 0;

  //(1, 1, 0) 1 in 1s
  cout << "numOfNNOr (matching 1, central): ";
  failed += compExp(numOfNNOr(5, 1), 6);

  //(1, 0, 0) 1 in 1s
  cout << "numOfNNOr (matching 1, boundary): ";
  failed += compExp(numOfNNOr(1, 1), 6);

  //(2, 1, 2) 2 in 2s
  cout << "numOfNNOr (matching 2, central): ";
  failed += compExp(numOfNNOr(38, 2), 6);

  //(2, 0, 2) 2 in 2s
  cout << "numOfNNOr (matching 2, boundary): ";
  failed += compExp(numOfNNOr(34, 2), 6);

  //(0, 3, 0) 2 in mixed
  cout << "numOfNNOr (perfect mix, boundary): ";
  failed += compExp(numOfNNOr(12, 2), 3);

  return failed;
}

int SimTest::are_nn() {
  int failed = 0;

  cout << "nn in i (central): ";
  failed += compExp(areNN(5, 6), 1);

  cout << "nn in j (boundary): ";
  failed += compExp(areNN(5, 9), 1);

  cout << "central, nn in k: ";
  failed += compExp(areNN(5, 53), 1);

  cout << "not nn: ";
  failed += compExp(areNN(0, 47), 0);

  cout << "not nn, but 1 away: ";
  failed += compExp(areNN(5, 7), 0);

  cout << "not nn, but diagonal: ";
  failed += compExp(areNN(5,10), 0);

  return failed;
}

//Test rotations where de = 0
int SimTest::rot_no_change() {
  int failed = 0;

  //(1, 1, 0)
  cout << "No change rotation (1 matching, central): ";
  failed += compExp(rotChange(5, 1), 0.0);

  //(1, 0, 0)
  cout << "No change rotation (1 matching, boundary): ";
  failed += compExp(rotChange(1, 1), 0.0);

  //(2, 2, 0)
  cout << "No change rotation (2 non-matching, central): ";
  failed += compExp(rotChange(10, 2), 0.0);

  //(3, 1, 0)
  cout << "No change rotation (2 non-matching, boundary): ";
  failed += compExp(rotChange(7, 2), 0.0);

  //(0, 3, 0)
  cout << "No change rotation (perfect mix, boundary): ";
  failed += compExp(rotChange(12, 2), 0.0);

  //(0, 3, 0) (3 orientations broken but 3 also made)
  cout << "No change rotation (perfect mix swap, boundary): ";
  failed += compExp(rotChange(12, 1), 0.0);

  return failed;
}

//Test rotations where de < 0
int SimTest::rot_good_change() {
  int failed = 0;

  //(1, 2, 2) 1 by all 2s -> 2 (6 orientations now match)
  cout << "Good change rotation (1->2, central): ";
  failed += compExp(rotChange(41, 2), -6.0);

  //(2, 2, 0) 2 by all 1s -> 1 (6 orientations now match)
  cout << "Good change rotation (2->1, central): ";
  failed += compExp(rotChange(10, 1), -6.0);

  //(3, 1, 0) 2 by all 1s -> 1 (6 orientations now match)
  cout << "Good change rotation (2->1, boundary): ";
  failed += compExp(rotChange(7, 1), -6.0);

  return failed;
}

//Test rotations where de > 0
int SimTest::rot_bad_change() {
  int failed = 0;

  //(1, 1, 0) 1 by all 1s -> 2 (6 orientations broken)
  cout << "Bad change rotation (1->2, central): ";
  failed += compExp(rotChange(5, 2), 6.0);

  //(1, 0, 0) 1 by all 1s -> 2 (6 orientations broken)
  cout << "Bad change rotation (1->2, boundary): ";
  failed += compExp(rotChange(1, 2), 6.0);

  //(2, 1, 2) 2 by all 2s -> 1 (6 orientations broken)
  cout << "Bad change rotation (2->1, central): ";
  failed += compExp(rotChange(38, 1), 6.0);

  //(2, 0, 2) 2 by all 2s -> 1 (6 orientations broken)
  cout << "Bad change rotation (2->1, boundary): ";
  failed += compExp(rotChange(34, 1), 6.0);

  return 0;
}

//Test swaps where de = 0
int SimTest::swap_no_change() {
  int failed = 0;

  //(3, 1, 0) with (2, 2, 0) 2s by all 1s, overlapping neighbors
  cout << "No change swap (non-matching 2s, overlapping neighbors, central/boundary): ";
  failed += compExp(swapChange(7, 10), 0.0);

  //(1, 1, 0) with (2, 2, 0) 1 by all 1s <-> 2 by all 1s (equal environments)
  cout << "No change swap (1 <-> 2 in all 1s, central): ";
  failed += compExp(swapChange(5, 10), 0.0);

  //(1, 0, 0) with (2, 2, 0) 1 by all 1s <-> 2 by all 1s (equal environments)
  cout << "No change swap (1 <-> 2 in all 1s, central/boundary): ";
  failed += compExp(swapChange(1, 10), 0.0);

  //(2, 2, 0) with (2, 1, 2) 2s by all 1s <-> 2s (equal orient, diff environments)
  cout << "No change swap (2s in all 1s <-> 2s, central): ";
  failed += compExp(swapChange(10, 38), 0.0);

  //(2, 2, 0) with (2, 0, 2) 2s by all 1s <-> 2s (equal orient)
  cout << "No change swap (2s in all 1s <-> 2s, central/boundary): ";
  failed += compExp(swapChange(10, 34), 0.0);

  return failed;
}

//Test swaps where de < 0
int SimTest::swap_good_change() {
  int failed = 0;

  //(1, 2, 2) with (2, 2, 0) 2 in all 1s <-> 1 in all 2s (max benefits (2*6*-1.5 = -18))
  cout << "Good change swap (2 in 1s <-> 1 in 2s, central): ";
  failed += compExp(swapChange(41, 10), -18.0);

  //(0, 2, 0) with (2, 1, 0) 2 in mostly 1s <-> 1 in mostly 2s (break 1+2, make 4+5 ((9-3)*-1.5 = -9)
  cout << "Good change swap (2 in mixed <-> 1 in 2s, central/boundary): ";
  failed += compExp(swapChange(8, 6), -9.0);

  return failed;
}

//Test swaps where de > 0
int SimTest::swap_bad_change() {
  int failed = 0;

  //(1, 1, 0) with (2, 1, 2) 1 in all 1s <-> 2 in all 2s (max breakage (2*6*1.5 = 18))
  cout << "Bad change swap (1 in 1s <-> 2 in 2s, central): ";
  failed += compExp(swapChange(5, 38), 18.0);

  return failed;
}

//Test neighbor swaps where de = 0
int SimTest::swap_nn_no_change() {
  int failed = 0;

  //(1, 1, 0) with (1, 0, 0) 1s surrounded by all 1s
  cout << "No change nn swap (matching 1s, central/boundary): ";
  failed += compExp(swapChange(5, 1), 0.0);

  //(2, 1, 2) with (2, 0, 2) 2s surrounded by all 2s
  cout << "No change nn swap (matching 2s, central/boundary): ";
  failed += compExp(swapChange(38, 34), 0.0);

  //(0, 0, 0) with (3, 0, 0) same identity, no change
  cout << "No change nn swap (1s, boundary): ";
  failed += compExp(swapChange(0, 3), 0.0);

  return failed;
}

//Test neighbor swaps where de < 0
int SimTest::swap_nn_good_change() {
  int failed = 0;

  //(1, 1, 2) with (1, 2, 2) 2 in mostly 1s <-> 1 in mostly 2s (break 2+0, make 5+3 ((8-2)*-1.5 = -9))
  cout << "Good change nn swap (2 in some 1s <-> 1 in some 2s, central): ";
  failed += compExp(swapChange(37, 41), -9.0);

  return failed;
}

//Test neighbor swaps where de > 0
int SimTest::swap_nn_bad_change() {
  int failed = 0;

  //(0, 3, 0) with (1, 3, 0) 2 in mixed <-> 1 in mostly 1s (break 3+5, make 0+2 ((2-8)*-1.5 = 9))
  cout << "Bad change nn swap (2 in mixed <-> 1 in some 1s, boundary): ";
  failed += compExp(swapChange(12, 13), 9);

  return failed;
}
