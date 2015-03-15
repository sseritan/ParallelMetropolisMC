//
//
// Written by Stefan Seritan on 03/06/15
// Header file for Simulation class
//

#ifndef _SIMULATION_HPP_
#define _SIMULATION_HPP_

//Parallel include
#include <tbb/flow_graph.h>

//Lightweight container to hold move info
class Move {
  public:
    Move(int ty, int p, int q) : type(ty), pos(p), par(q) {};
    int type; //0 for rotation, 1 for particle swap
    int pos; //Lattice position
    int par; //Parameter: will be new or for rotation, 2nd position if swap
};

//Simulation class
//Holds lattice, generates and performs moves, calculates data
class Simulation {
    //Hardcoded simulation parameters
    static const int Q = 6;
    static const double K = 0.5;
    static const double A = 1.0;
    static const double ROTATION = 0.5;
    static const double PARTSWAP = 0.5;

    //Inputted simulation parameters
    const int Lx, Ly, Lz; //3D Lattice dimensions
    const int NMAX; //Shortcut for Lx*Ly*Lz
    const double kT; //Temperature
    const double cutoff; //Theta cutoff value for phase separation

    double energy; //Energy (continuously updated)

    int* array; //1D array of id and orientation
    //tens place is id, ones place is orientation

    //Private functions
    //Energy functions
    double rotChange(int pos, int q);
    double swapChange(int pos1, int pos2);
    //Move functions
    void performMove(Move* m);
    //Data functions
    double* calctheta();
    double* calcTheta();
    //Cell query functions
    double pairEnergy(int pos1, int pos2);
    int numOfNNId(int pos, int i);
    int numOfNNOr(int pos, int o);
    int areNN(int pos1, int pos2);
    //Utility functions
    int wrap1d(int coord, int dir, int step);
    int step3d(int index, int dir, int step);

    //Helper struct for dependency graph
    struct moveBody {
      Simulation* sim;
      Move* move;
      moveBody(Simulation* s, Move* m) : sim(s), move(m) {};
      void operator() (tbb::flow::continue_msg) const {
        sim->performMove(move);
      }
    };

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
