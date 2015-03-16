//
// Written by Stefan Seritan on 03/11/15 for CS 140
// Modified by Wei Dai 03/15/15
// Cell and Move classes (encapsulation for Simulation data)
//

#ifndef _CELL_MOVE_HPP_
#define _CELL_MOVE_HPP_

#include <iostream>
#include <atomic>

class SpinLock {
    std::atomic_flag flag;
  public:
    void lock() {
      while(flag.test_and_set(std::memory_order_acquire));
    }
    void unlock() {
      flag.clear(std::memory_order_release);
    }
};

//Cell class
//Store id, orientation, and update history
class Cell: public SpinLock {
    int id, orient; //Identity and orientation

  public:
    Cell(){};
    Cell(int i, int o) : id(i), orient(o) {};
    //Getters
    int getId() const {return id;}
    int getOr() const {return orient;}
    //Setters
    void setId(int i) {id = i;}
    void setOr(int o) {orient = o;}
    //Utility functions
    void printCell() const {std::cout << "Id " << id << " Or " << orient << std::endl;}
};

//Lightweight container to hold move info
class Move {
    const int type; //0 for rotation, 1 for particle swap
    const int pos; //Lattice position
    const int par; //Parameter: will be new or for rotation, 2nd position if swap
    const double prob; //a random number in [0, 1]

  public:
    Move(int ty, int p, int q, double r) : type(ty), pos(p), par(q), prob(r) {};
    //Default destructor is ok
    //Getters
    int getType() const {return type;}
    int getPos() const {return pos;}
    int getPar() const {return par;}
    double getProb() const {return prob;}
    void printMove() const {
      if (type == 0) {
        std::cout << "Rot(" << pos << ", " << par << ")" << " with prob = "  << prob << std::endl;
      } else if (type == 1) {
        std::cout << "Swap(" << pos << ", " << par << ")" << " with prob = "  << prob << std::endl;
      }
    }
    //No setters
};

#endif
