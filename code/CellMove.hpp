//
// Written by Stefan Seritan on 03/11/15 for CS 140
// Cell and Move classes (encapsulation for Simulation data)
//

#ifndef _CELL_MOVE_HPP_
#define _CELL_MOVE_HPP_

#include <iostream>

//Cell class
//Store id, orientation, and update history
class Cell {
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
    int type; //0 for rotation, 1 for particle swap
    int pos; //Lattice position
    int par; //Parameter: will be new or for rotation, 2nd position if swap

  public:
    Move(int ty, int p, int q) : type(ty), pos(p), par(q) {};
    //Default destructor is ok
    //Getters
    int getType() const {return type;}
    int getPos() const {return pos;}
    int getPar() const {return par;}
    //No setters
};

#endif
