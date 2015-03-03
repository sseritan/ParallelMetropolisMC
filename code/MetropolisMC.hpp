//
//
// Originally written by Stefan Seritan on 8/28/14 for the Peters Group
// Modified by Stefan Seritan and Wei Dai for CS 140 Final Project Winter 2015
//
// Header file for MetropolisMC and subroutine
// Includes simulation parameters, templated data structures, and function prototypes
//
//

#ifndef _METROPOLIS_MC_HPP_
#define _METROPOLIS_MC_HPP_

//Array include
#include <array>
//Vector include
#include <vector>

/* SIMULATION PARAMETERS */

//Orientations
#define Q 6

//Lattice Size
#define Lx 10
#define Ly 10
#define Lz 10

//Hamiltonian definition
#define K 0.5
#define A 1.0

//Move probabilities
#define ROTATION 0.5
#define PARTSWAP 0.5

#define np 4

/* TEMPLATED DATA STRUCTURES */

// datatype representing a move
struct Move {
  char type; // 0 for rotation, 1 for swap
  int cell0[3];
  int cell1[3];
};

//Template class to make a 3d array
template <typename T, size_t x, size_t y, size_t z>
using Dim3Array = std::array<std::array<std::array<T, z>, y>, x>;

//Specific 3D array that is the fixed size of the box
template <typename T>
using SimArray = Dim3Array<T, Lx, Ly, Lz>;

//Specific 2D array for the 2d histogram
using Dim2Array = std::array<std::array<double, 100>, 100>;

/* FUNCTION PROTOTYPES */

void gen_moves(SimArray<int>& Y, std::vector<Move>* M, int& moves_left);
int check_conflicts(SimArray<int>& Y, Move* move, int pc);
void print_move(Move& Move);

//Main simulation functions
void sweep(std::vector<Move>* M, SimArray<int>& X, SimArray<int>& S, double& e, const double kT);
int mod(int k, int n);

//Energy related functions
double energy(const SimArray<int>& X, const SimArray<int>& S);
double pairwise_energy(int m1, int m2, int s1, int s2);
double point_energy(const SimArray<int>& X, const SimArray<int>& S, int i, int j, int k);
double rotation_energy_change(const SimArray<int>& S, int i, int j, int k, int q2);
double particle_swap_energy_change(SimArray<int>& X, SimArray<int>& S, int i, int j, int k, int ii, int jj, int kk);

//Order parameter related functions
SimArray<double> phase_parameter(const SimArray<int>& X);
double local_phase_parameter(const SimArray<int>& X, int i, int j, int k);
double medium_phase_parameter(const SimArray<double>& theta, int i, int j, int k);

//Data collection related functions
std::array<double, 2> phase_data(const SimArray<int>& X, const SimArray<double>& Theta, const double cutoff);

//I/O Functions
void print_VMD_snapshot(SimArray<int>& X, SimArray<int>& S, int t);

#endif

