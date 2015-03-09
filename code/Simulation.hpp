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
    Cell(int i, int o);
    //Getters
    int getId() {return id;}
    int getOr() {return orient;}
    int getLastUpdate() {return updateHistory.back();}
    //Setters
    void setId(int i) {id = i;}
    void setOr(int o) {orient = o;}
    void pushUpdate(int m);
    void resetUpdate();
    //Utility functions
    void printCell();
};

//Simulation
class Simulation {
    int Lx, Ly, Lz; //3D Lattice dimensions
    int NMAX;
    double kT; //Temperature
    double energy; //Energy (continuously updated)
    double cutoff; //Theta cutoff value for phase separation

    Cell** array; //1D array of Cell*
    //std::atomic_flag* flags; //Array of flags (have atomic test and set)

    //Private functions
    //Energy functions
    double rotChange(int pos, int q);
    double swapChange(int pos1, int pos2);
    //Data functions
    double* calctheta();
    double* calcTheta();
    //Periodic boundary functions
    int wrap1d(int coord, int dir, int step);
    //Cell query functions
    double pairEnergy(int pos1, int pos2);
    int numOfNNId(int pos, int i);
    int numOfNNOr(int pos, int o);
    int areNN(int pos1, int pos2);
    //Cell manipulation functions
    void swapIdOr(Cell* c1, Cell* c2);
    //Utility functions
    int step3d(int index, int dir, int step);

    //Friend for testing
    friend class SimTest;
  public:
    //Constructor
    Simulation(int x, int y, int z, double T, double compA, double c);
    //Destructor
    ~Simulation();
    //Move functions
    void doSweep();
    //Data functions
    double getEnergy();
    double* calcThetaHistogram();
    double* calcX1();
};

#endif
