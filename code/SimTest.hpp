//
// Written by Stefan Seritan on 03/08/15 for CS 140
// Header file for SimTest class
//

#ifndef _SIM_TEST_HPP_
#define _SIM_TEST_HPP_

#include "./Simulation.hpp"

class SimTest : public Simulation {
  public:
    SimTest();
    //Default destructor is enough
    //TESTING FUNCTIONS
    //Test helper function
    int testChange(double r, double e);
    //Rotation tests
    int rot_no_change();
    int rot_good_change();
    int rot_bad_change();
    //Swap tests
    int swap_no_change();
    int swap_good_change();
    int swap_bad_change();
    int swap_nn_no_change();
    int swap_nn_good_change();
    int swap_nn_bad_change();
};

#endif
