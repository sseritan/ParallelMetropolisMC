//
// Written by Stefan Seritan on 03/08/15
// Test suite for Simulation rotChange and dataChange
//

#include <iostream>

#include "./CellMove.hpp"
#include "./Simulation.hpp"
#include "./SimTest.hpp"

using namespace std;

int main() {
  //Initialize my own special mini simulation (5x5x5)
  SimTest* st = new SimTest();

  int failed = 0;

  //Test periodic boundary functions
  cout << "PERIODIC BOUNDARY TESTS:" << endl;
  failed += st->wrap();
  failed += st->step();

  //Test neighbor functions
  cout << "NEIGHBOR TESTS:" << endl;
  failed += st->nn_id_count();
  failed += st->nn_or_count();
  failed += st->are_nn();

  //Test rotChange
  cout << "ROTATION TESTS:" << endl;
  failed += st->rot_no_change();
  failed += st->rot_good_change();
  failed += st->rot_bad_change();

  //Test swapChange
  cout << "SWAP TESTS:" << endl;
  failed += st->swap_no_change();
  failed += st->swap_good_change();
  failed +=st->swap_bad_change();

  //Test neighbor swap in swapChange (tricky cases)
  cout << "SWAP NEIGHBORS TESTS:" << endl;
  failed += st->swap_nn_no_change();
  failed += st->swap_nn_good_change();
  failed += st->swap_nn_bad_change();

  if (failed) {
    cout << "\nFailed " << failed << " tests. :(" << endl;
  } else {
    cout << "\nPassed all tests! :D Here is a congratulatory sweet potato 🍠" << endl;
  }

  //Memory Management
  delete st;
}
