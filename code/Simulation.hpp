//
//
// Written by Stefan Seritan on 03/06/15
// Modified by Wei Dai 03/15/15
// Header file for Simulation class
//

#ifndef _SIMULATION_HPP_
#define _SIMULATION_HPP_

#include "./CellMove.hpp"
#include <random>

#define CACHE_LINE_SIZE 64

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
    mutable double energy; //Energy (continuously updated)
    double cutoff; //Theta cutoff value for phase separation

    mutable Cell *array; //1D array of Cell*
    /* alignas(CACHE_LINE_SIZE) mutable Cell *array; //1D array of Cell* */

    //Private functions
    //Locking
    void posLocks(int pos, int& even, int& odd) const;
    //Energy functions
    double rotChange(int pos, int q) const;
    double swapChange(int pos1, int pos2) const;
    double performMove(const Move* const m) const;
    //Data functions
    double* calctheta() const;
    double* calcTheta() const;
    //Periodic boundary functions
    int wrap1d(int coord, int dir, int step) const;
    //Cell query functions
    double pairEnergy(int pos1, int pos2) const;
    int numOfNNId(int pos, int i) const;
    int numOfNNOr(int pos, int o) const;
    int areNN(int pos1, int pos2) const;
    //Cell manipulation functions
    void swapIdOr(Cell& c1, Cell& c2) const;
    //Utility functions
    int step3d(int index, int dir, int step) const;

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
    double getEnergy() const;
    double* calcThetaHistogram() const;
    double* calcX1() const;
};

#endif
