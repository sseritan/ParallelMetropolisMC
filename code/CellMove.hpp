//
// Written by Stefan Seritan on 03/11/15 for CS 140
// Cell and Move classes (encapsulation for Simulation data)
//

#ifndef _CELL_MOVE_HPP_
#define _CELL_MOVE_HPP_

//Cell class
//Store id, orientation, and update history
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

//Lightweight container to hold move info
class Move {
    int type; //0 for rotation, 1 for particle swap
    int pos; //Lattice position
    int par; //Parameter: will be new or for rotation, 2nd position if swap

  public:
    Move(int ty, int p, int q) : type(ty), pos(p), par(q) {};
    //Default destructor is ok
    //Getters
    int getType() {return type;}
    int getPos() {return pos;}
    int getPar() {return par;}
};

#endif
