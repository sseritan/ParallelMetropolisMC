//
//
// Originally written by Stefan Seritan on 8/28/14 for the Peters Group
// Modified by Stefan Seritan and Wei Dai for CS 140 Final Project Winter 2015
//
// Header file for MetropolisMC and subroutine
// Includes simulation parameters and templated data structures
//
//

#ifndef _METROPOLIS_MC_HPP_
#define _METROPOLIS_MC_HPP_

//Array include
#include <array>

/* SIMULATION PARAMETERS */

//Temperature
#define kT 3.0

//Orientations
#define Q 6

//Hamiltonian definition
#define K11 0.5
#define K22 0.4
#define A 3.0

//Box size
#define Lx 10
#define Ly 10
#define Lz 10

//Move probabilities
#define ROTATION 0.5
#define PARTSWAP 0.5

//Composition (total fraction that is A)
#define COMPA 1.0

//Sweep information
#define EQ_SWEEP 5000
#define DATA_SWEEP 10000

//Cutoff value for Theta to determine phases in a Solid-Solid Phase Diagram Run
#define THETA_CUTOFF 0.5

//Cutoff value for Phi to determine phases in a Liquid-Solid Phase Diagram Run
#define PHI_CUTOFF 0.5

//Parameter that determines what info is kept
//1: MeltingPoint run (keep energy, COMP should be 0.0 or 1.0)
//2: Cutoff run (keep 2d histogram to find cutoff values that divide phases)
//3: Solid-Solid Phase Diagram run (use THETA_CUTOFF to get XA for diff. phases)
//4: Liquid-Solid Phase Diagram run (use PHI_CUTOFF to get XA for diff. phases)
#define RUNTYPE 2

//Debugging flag will output a lot of information about individual moves
//Recommend have only one sweep with a small box
#define DEBUGGING 0

/* TEMPLATED DATA STRUCTURES */

//Template class to make a 3d array
template <typename T, size_t x, size_t y, size_t z>
using Dim3Array = std::array<std::array<std::array<T, z>, y>, x>;

//Specific 3D array that is the fixed size of the box
template <typename T>
using SimArray = Dim3Array<T, Lx, Ly, Lz>;

//Specific 2D array for the 2d histogram
using Dim2Array = std::array<std::array<float, 100>, 100>;

#endif

