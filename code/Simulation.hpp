//
//
// Written by Stefan Seritan on 03/06/15
// Header file for Simulation class
//

#ifndef _SIMULATION_HPP_
#define _SIMULATION_HPP_

#include "./CellMove.hpp"
#include <mutex>

//Simulation class
//Holds lattice, generates and performs moves, calculates data
class Simulation {
    //Hardcoded simulation parameters
    const int Q = 6;
    const double K = 0.5;
    const double A = 1.0;
    const double ROTATION = 0.5;
    const double PARTSWAP = 0.5;

    int Lx, Ly, Lz; //3D Lattice dimensions
    int NMAX;
    double kT; //Temperature
    double energy; //Energy (continuously updated)
    double cutoff; //Theta cutoff value for phase separation

    Cell** array; //1D array of Cell*

    //std::atomic_flag* locks; // Even and odd locks (have atomic test and set)
    std::mutex* locks; // Even and odd locks from STL

    //Private functions
    //Locking
    void posLocks(int pos, int& even, int& odd);
    //Energy functions
    double rotChange(int pos, int q);
    double swapChange(int pos1, int pos2);
    void performMove(Move* m);
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
