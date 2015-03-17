#!/bin/bash

#Written by Stefan Seritan on 03/02/15 for CS 140
#Quick script that runs the MetropolisMC simulation to get melting temp information
#Melting temp is given by peak in heat capacity, which is dE/dt

KT=1.00
STEPSIZE=0.05
MAX=1.50

STEPS=$(echo "(($MAX - $KT)/$STEPSIZE)+1" | bc)
COUNTER=0

while [ $COUNTER -lt $STEPS ]; do
  echo "Running "$KT

  #Run MC simulation at this temperature
  FILENAME=$(echo "meltingTemp_"$KT".txt")
  mpirun -machinefile $PBS_NODEFILE -cpus-per-proc 4 -npernode 2 ./MC $KT 10 10 100 5000 10000 1.0 1.0 > $FILENAME

  #Set up next iteration
  KT=$(echo "x=$KT+$STEPSIZE; if(x<1) print 0; x" | bc)
  let COUNTER=COUNTER+1
done
