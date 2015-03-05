#!/bin/bash

#Written by Stefan Seritan on 03/02/15 for CS 140
#Quick script that runs the MetropolisMC simulation to get timing information

#Run increasing sweep number on constant procs (16)
for i in 100 500 1000 5000 10000; do
  echo "Running "$i" on 16 procs"

  #Run MC simulation at this temperature
  FILENAME=$(echo "scaling_"$i".txt")
  CILK_NWORKERS=16 ./MetropolisMC 1.0 $i $i 0.5 0.5 > $FILENAME
done

#Run constant sweep size (1000) on increasing procs
for p in 1 2 4 6 8 12 16; do
  echo "Running 5000 on "$p" procs"

  FILENAME=$(echo "proc_"$p".txt")
  CILK_NWORKERS=$p ./MetropolisMC 1.0 5000 5000 0.5 0.5 > $FILENAME
done
