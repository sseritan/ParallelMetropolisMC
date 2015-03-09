//
// Written by Stefan Seritan on 03/08/15 for CS 140
// Contains the SimTest class, used for testing Simulation class
//

#include <iostream>

#include "./MetropolisMC.hpp"
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
  for (int k = 0; k < 4; k++) {
    for (int j = 0; j < 4; j++) {
      for (int i = 0; i < 4; i++) {
        cout << array[i+ Lx*j + Lx*Ly*k]->getId() << array[i+ Lx*j + Lx*Ly*k]->getOr() << " ";
      }
      cout << endl;
    }
    cout << endl;
  }

  //Initialize energies
  energy = 0.0;
  for (int i = 0; i < NMAX; i++) {
    //Calculate pairwise energy with forward pairs in each direction
    for (int j = 0; j < 3; j++) {
      energy += pairEnergy(i, step3d(i, j, 1));
    }
  }
}

int SimTest::testChange(double r, double e) {
  if (r != e) {
    cout << "Failed. Expected " << e << " but got " << r << "." << endl;
    return 1;
  }
  cout << "Passed!" << endl;
  return 0;
}

//Test rotations where de = 0
int SimTest::rot_no_change() {
  int failed = 0;

  //(1, 1, 0)
  cout << "No change rotation (1 matching, central): ";
  failed += testChange(rotChange(5, 1), 0.0);

  //(1, 0, 0)
  cout << "No change rotation (1 matching, boundary): ";
  failed += testChange(rotChange(1, 1), 0.0);

  //(2, 2, 0)
  cout << "No change rotation (2 non-matching, central): ";
  failed += testChange(rotChange(10, 2), 0.0);

  //(3, 1, 0)
  cout << "No change rotation (2 non-matching, boundary): ";
  failed += testChange(rotChange(7, 2), 0.0);

  //(0, 3, 0)
  cout << "No change rotation (perfect mix, boundary): ";
  failed += testChange(rotChange(12, 2), 0.0);

  //(0, 3, 0) (3 orientations broken but 3 also made)
  cout << "No change rotation (perfect mix swap, boundary): ";
  failed += testChange(rotChange(12, 1), 0.0);

  return failed;
}

//Test rotations where de < 0
int SimTest::rot_good_change() {
  int failed = 0;

  //(1, 2, 2) 1 by all 2s -> 2 (6 orientations now match)
  cout << "Good change rotation (1->2, central): ";
  failed += testChange(rotChange(41, 2), -6.0);

  //(2, 2, 0) 2 by all 1s -> 1 (6 orientations now match)
  cout << "Good change rotation (2->1, central): ";
  failed += testChange(rotChange(10, 1), -6.0);

  //(3, 1, 0) 2 by all 1s -> 1 (6 orientations now match)
  cout << "Good change rotation (2->1, boundary): ";
  failed += testChange(rotChange(7, 1), -6.0);

  return failed;
}

//Test rotations where de > 0
int SimTest::rot_bad_change() {
  int failed = 0;

  //(1, 1, 0) 1 by all 1s -> 2 (6 orientations broken)
  cout << "Bad change rotation (1->2, central): ";
  failed += testChange(rotChange(5, 2), 6.0);

  //(1, 0, 0) 1 by all 1s -> 2 (6 orientations broken)
  cout << "Bad change rotation (1->2, boundary): ";
  failed += testChange(rotChange(1, 2), 6.0);

  //(2, 1, 2) 2 by all 2s -> 1 (6 orientations broken)
  cout << "Bad change rotation (2->1, central): ";
  failed += testChange(rotChange(38, 1), 6.0);

  //(2, 0, 2) 2 by all 2s -> 1 (6 orientations broken)
  cout << "Bad change rotation (2->1, boundary): ";
  failed += testChange(rotChange(34, 1), 6.0);

  return 0;
}

//Test swaps where de = 0
int SimTest::swap_no_change() {
  int failed = 0;

  //(3, 1, 0) with (2, 2, 0) 2s by all 1s, overlapping neighbors
  cout << "No change swap (non-matching 2s, overlapping neighbors, central/boundary): ";
  failed += testChange(swapChange(7, 10), 0.0);

  //(1, 1, 0) with (2, 2, 0) 1 by all 1s <-> 2 by all 1s (equal environments)
  cout << "No change swap (1 <-> 2 in all 1s, central): ";
  failed += testChange(swapChange(5, 10), 0.0);

  //(1, 0, 0) with (2, 2, 0) 1 by all 1s <-> 2 by all 1s (equal environments)
  cout << "No change swap (1 <-> 2 in all 1s, central/boundary): ";
  failed += testChange(swapChange(1, 10), 0.0);

  //(2, 2, 0) with (2, 1, 2) 2s by all 1s <-> 2s (equal orient, diff environments)
  cout << "No change swap (2s in all 1s <-> 2s, central): ";
  failed += testChange(swapChange(10, 38), 0.0);

  //(2, 2, 0) with (2, 0, 2) 2s by all 1s <-> 2s (equal orient)
  cout << "No change swap (2s in all 1s <-> 2s, central/boundary): ";
  failed += testChange(swapChange(10, 34), 0.0);

  return failed;
}

//Test swaps where de < 0
int SimTest::swap_good_change() {
  int failed = 0;

  //(1, 2, 2) with (2, 2, 0) 2 in all 1s <-> 1 in all 2s (max benefits (2*6*-1.5 = -18))
  cout << "Good change swap (2 in 1s <-> 1 in 2s, central): ";
  failed += testChange(swapChange(41, 10), -18.0);

  //(0, 2, 0) with (2, 1, 0) 2 in mostly 1s <-> 1 in mostly 2s (break 1+2, make 4+5 ((9-3)*-1.5 = -9)
  cout << "Good change swap (2 in mixed <-> 1 in 2s, central/boundary): ";
  failed += testChange(swapChange(8, 6), -9.0);

  return failed;
}

//Test swaps where de > 0
int SimTest::swap_bad_change() {
  int failed = 0;

  //(1, 1, 0) with (2, 1, 2) 1 in all 1s <-> 2 in all 2s (max breakage (2*6*1.5 = 18))
  cout << "Bad change swap (1 in 1s <-> 2 in 2s, central): ";
  failed += testChange(swapChange(5, 38), 18.0);

  return failed;
}

//Test neighbor swaps where de = 0
int SimTest::swap_nn_no_change() {
  int failed = 0;

  //(1, 1, 0) with (1, 0, 0) 1s surrounded by all 1s
  cout << "No change nn swap (matching 1s, central/boundary): ";
  failed += testChange(swapChange(5, 1), 0.0);

  //(2, 1, 2) with (2, 0, 2) 2s surrounded by all 2s
  cout << "No change nn swap (matching 2s, central/boundary): ";
  failed += testChange(swapChange(38, 34), 0.0);

  //(0, 0, 0) with (3, 0, 0) same identity, no change
  cout << "No change nn swap (1s, boundary): ";
  failed += testChange(swapChange(0, 3), 0.0);

  return failed;
}

//Test neighbor swaps where de < 0
int SimTest::swap_nn_good_change() {
  int failed = 0;

  //(1, 1, 2) with (1, 2, 2) 2 in mostly 1s <-> 1 in mostly 2s (break 2+0, make 5+3 ((8-2)*-1.5 = -9))
  cout << "Good change nn swap (2 in some 1s <-> 1 in some 2s, central): ";
  failed += testChange(swapChange(37, 41), -9.0);

  return failed;
}

//Test neighbor swaps where de > 0
int SimTest::swap_nn_bad_change() {
  int failed = 0;

  //(0, 3, 0) with (1, 3, 0) 2 in mixed <-> 1 in mostly 1s (break 3+5, make 0+3 ((3-8)*-1.5 = 7.5))
  cout << "Bad change nn swap (2 in mixed <-> 1 in some 1s, boundary): ";
  failed += testChange(swapChange(12, 13), 7.5);

  return failed;
}
