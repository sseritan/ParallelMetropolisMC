//
// Written by Stefan Seritan on 03/06/15 for CS 140
// Simulation.cpp
// All member functions for the Cell and Simulation classes
//

#include <iostream>
#include <vector>

//Local include
#include "Simulation.hpp"

using namespace std;

/* Cell FUNCTIONS */
//Getters
int Cell::getId() {
  return id;
}

int Cell::getOrientation() {
  return orient;
}

int Cell::getLastUpdate() {
  return updateHistory.back();
}

//Setters
void Cell::updateCell(int i, int o, int m) {
  id = i;
  orient = o;
  updateHistory.push_back(m);
}

void Cell::resetUpdate() {
  updateHistory.clear();
}

/* Simulation FUNCTIONS */
//Constructor
Simulation::Simulation(int x, int y, int z, double compA) {
  Lx = x; Ly = y; Lz = z;

  //Allocate array of Cells
  cout << "Allocated simulation array of " << Lx << "x" << Ly << "x" << Lz << endl;
  array = new Cell [Lx*Ly*Lz];

  //Initialize cells
  int threshold = compA*Lx*Ly*Lz;
  cout << "Simulation seeded with " << threshold << " particles of species 1" << endl;

  //Set species 1
  for (int i = 0; i < threshold; i++) {
    array[i].updateCell(1, 1, 0);
  }
  //Set species 2
  for (int i = threshold; i < Lx*Ly*Lz; i++) {
    array[i].updateCell(2, 1, 0);
  }
}

//Destructor
Simulation::~Simulation() {
  //Memory Management
  delete[] array;
}

//Operator () overloading
Cell* Simulation::operator() (int x, int y, int z) {
  return &array[x+Lx*y+Lx*Ly*z];
}
