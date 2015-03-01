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

/* SIMULATION PARAMETERS */

//Orientations
#define Q 6

//Lattice Size
#define Lx 10
#define Ly 10
#define Lz 100

//Hamiltonian definition
#define K 0.5
#define A 1.0

//Move probabilities
#define ROTATION 0.5
#define PARTSWAP 0.5

/* TEMPLATED DATA STRUCTURES */

//Template class to make a 3d array
template <typename T, size_t x, size_t y, size_t z>
using Dim3Array = std::array<std::array<std::array<T, z>, y>, x>;

//Specific 3D array that is the fixed size of the box
template <typename T>
using SimArray = Dim3Array<T, Lx, Ly, Lz>;

//Specific 2D array for the 2d histogram
using Dim2Array = std::array<std::array<float, 100>, 100>;

/* FUNCTION PROTOTYPES */
//Main simulation functions
void sweep(SimArray<int>& X, SimArray<int>& S, double& e, const double kT);
int mod(int k, int n);

//Energy related functions
double energy(const SimArray<int>& X, const SimArray<int>& S);
double pairwise_energy(int m1, int m2, int s1, int s2);
double point_energy(const SimArray<int>& X, const SimArray<int>& S, int i, int j, int k);
double rotation_energy_change(const SimArray<int>& X, SimArray<int>& S, int i, int j, int k, int q);
double particle_swap_energy_change(SimArray<int>& X, SimArray<int>& S, int i, int j, int k, int ii, int jj, int kk);

//Order parameter related functions
SimArray<float> phase_parameter(const SimArray<int>& X);
float local_phase_parameter(const SimArray<int>& X, int i, int j, int k);
float medium_phase_parameter(const SimArray<float>& theta, int i, int j, int k);
SimArray<float> orientation_parameter(const SimArray<int>& S);
float local_orientation_parameter(const SimArray<int>& S, int i, int j, int k);
float medium_orientation_parameter(const SimArray<float>& phi, int i, int j, int k);

//Data collection related functions
std::array<float, 100> histogram(const SimArray<float>& param);
Dim2Array histogram2d(const SimArray<float>& Theta, const SimArray<float>& Phi);
std::array<float, 2> phase_data(const SimArray<int>& X, const SimArray<float>& param, const double cutoff);


//I/O Functions
int parseArgs(int& rT, double& kT, int& eS, int& dS, double& cA, double& pC);
void print_VMD_snapshot(SimArray<int>& X, SimArray<int>& S, int t);

#endif

