#!/bin/bash

#Written by Stefan Seritan on 03/02/15 for CS 140
#Quick script that runs the MetropolisMC simulation to get timing information

#Run increasing sweep number on constant procs (1)
for i in 100 500 1000 5000 10000; do
  echo "Running "$i" on 1 proc (because serial)"

  #Run MC simulation at this temperature
  FILENAME=$(echo "scaling_"$i".txt")
  ./MC 1.0 10 10 100 $i $i 0.5 0.5 > $FILENAME
done
