//
// Written by Stefan Seritan on 03/08/15
// Test suite for Simulation rotChange and dataChange
//

#include <iostream>

//Testing hack to have access to change functions
#define private public

#include "./MetropolisMC.hpp"
#include "./Simulation.hpp"

using namespace std;

int main() {
  //Initialize my own special mini simulation (5x5x5)
  Simulation* sim = new Simulation(-1);

  int failed = 0;

  //Test rotChange
  cout << "ROTATION TESTS:" << endl;
  failed += sim->rot_no_change();
  failed += sim->rot_good_change();
  failed += sim->rot_bad_change();

  //Test swapChange
  cout << "SWAP TESTS:" << endl;
  failed += sim->swap_no_change();
  failed += sim->swap_good_change();
  failed +=sim->swap_bad_change();

  //Test neighbor swap in swapChange (tricky cases)
  cout << "SWAP NEIGHBORS TESTS:" << endl;
  failed += sim->swap_nn_no_change();
  failed += sim->swap_nn_good_change();
  failed += sim->swap_nn_bad_change();

  if (failed) {
    cout << "\nFailed " << failed << " tests. :(" << endl;
  } else {
    cout << "\nPassed all tests! :D" << endl;
  }

  //Memory Management
  delete sim;
}
