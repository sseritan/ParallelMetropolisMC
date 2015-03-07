//
//
// Written by Stefan Seritan on 03/06/15
// Header file for Simulation.cpp
//

#ifndef _SIMULATION_HPP_
#define _SIMULATION_HPP_

#include <vector>
#include <atomic>

//Cell class
//Store id, orientation, and update history
class Cell {
    int id, orient; //Identity and orientation
    std::vector<int> updateHistory; //Array of ints to keep track of which moves update the cell
  public:
    //Getters
    int getId();
    int getOrientation();
    int getLastUpdate();
    //Setters
    void updateCell(int i, int o, int m);
    void resetUpdate();
};

//Simulation
class Simulation {
    int Lx, Ly, Lz;
    Cell* array; //1D array of Cells
    //std::atomic_flag* flags; //Array of flags (have atomic test and set)
  public:
    //Constructor
    Simulation(int x, int y, int z, double compA);
    //Destructor
    ~Simulation();
    //Operator (x,y,z) overload
    Cell* operator() (int x, int y, int z);
};

#endif
