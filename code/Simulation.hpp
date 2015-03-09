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
    Cell* im; Cell* ip; //Neighbors in the i direction
    Cell* jm; Cell* jp; //Neighbors in the j direction
    Cell* km; Cell* kp; //Neighbors in the k direction
    double theta; //Local value of theta

  public:
    Cell(int i, int o);
    //Getters
    int getId() {return id;}
    int getOr() {return orient;}
    int getLastUpdate() {return updateHistory.back();}
    Cell* getIm() {return im;}
    Cell* getIp() {return ip;}
    Cell* getJm() {return jm;}
    Cell* getJp() {return jp;}
    Cell* getKm() {return km;}
    Cell* getKp() {return kp;}
    //Setters
    void setId(int i) {id = i;}
    void setOr(int o) {orient = o;}
    void pushUpdate(int m);
    void resetUpdate();
    void setIm(Cell *c) {im = c;}
    void setIp(Cell *c) {ip = c;}
    void setJm(Cell *c) {jm = c;}
    void setJp(Cell *c) {jp = c;}
    void setKm(Cell *c) {km = c;}
    void setKp(Cell *c) {kp = c;}
    //Energy functions
    double pairEnergy(int i, int o, int q);
    double pointEnergy(int i, int o);
    //Data functions
    void thetaInit();
    double calcTheta();
    //Utility functions
    int isNeighbor(Cell* c);
    void swapIdOr(Cell* c);
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
    int wrap1d(int coord, int dir, int step);
    int wrap3d(int index, int dir, int step);
    double rotChange(Cell* c, int q);
    double swapChange(Cell* c1, Cell* c2);
    double* calctheta();
    double* calcTheta();

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
