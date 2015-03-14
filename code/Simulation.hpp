//
//
// Written by Stefan Seritan on 03/06/15
// Header file for Simulation class
//

#ifndef _SIMULATION_HPP_
#define _SIMULATION_HPP_

//Cell class
//Store id, orientation, and update history
//TODO: Move away from Cell?
class Cell {
    int id, orient; //Identity and orientation

  public:
    Cell(int i, int o) : id(i), orient(o) {};
    //Getters
    int getId() {return id;}
    int getOr() {return orient;}
    //Setters
    void setId(int i) {id = i;}
    void setOr(int o) {orient = o;}
    //Utility functions
    void printCell() {std::cout << "Id " << id << " Or " << orient << std::endl;}
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

    Cell** array; //1D array of Cell*

    //Private functions
    //Energy functions
    double rotChange(int pos, int q);
    double swapChange(int pos1, int pos2);
    //Move functions
    void performMove(int type, int pos, int par);
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
