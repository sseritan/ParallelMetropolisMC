//
// Written by Stefan Seritan on 03/10/15 for CS 140
// Header file for Move and MoveQueue
//

#ifndef _MOVE_THREAD_HPP_
#define _MOVE_THREAD_HPP_

//Lightweight container to hold move info
class Move {
    int type; //0 for rotation, 1 for particle swap
    int time; //time id, used for conflict resolution
    int pos; //Lattice position
    int par; //Parameter: will be new or for rotation, 2nd position if swap

  public:
    Move(int ty, int ti, int p, int q);
    //Default destructor is ok
    //Getters
    int getType() {return type;}
    int getTime() {return time;}
    int getPos() {return pos;}
    int getPar() {return par;}
    //No setters
};

#endif
