# ParallelMetropolisMC

The original version of this code was written by Stefan Seritan on
August 28th, 2014 for his work on solid-solid equilibria for the Peters
Group in the Department of Chemical Engineering at UCSB.

The code was used as a basic template of a Metropolis Monte Carlo
simulation for a Potts Lattice Gas model of a binary mixture. For CS
140 (Parallel Scientific Computing), the code was parallelized as a
final project.

Basic Goals for the CS 140 Final Project:
- First streamline the Metropolis Monte Carlo simulation to make a more
  efficient algorithm overall
- Successfully parallelize the Monte Carlo simulation while maintaing
  the correctness of the simulation results, as compared with the
  original serial version.
- If possible, test several ideas for parallelization and report on
  performance details.
