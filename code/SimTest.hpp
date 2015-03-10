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
    //Periodic boundary tests
    int wrap();
    int step();
    //Neighbor info tests
    int nn_id_count();
    int nn_or_count();
    int are_nn();
    //Test helper function
    int compExp(double r, double e);
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
