#!/bin/bash

#Written by Stefan Seritan on 03/02/15 for CS 140
#Quick script that runs the MetropolisMC simulation to get melting temp information
#Melting temp is given by peak in heat capacity, which is dE/dt

KT=0.5
STEPSIZE=0.05
MAX=2.0

STEPS=$(echo "(($MAX - $KT)/$STEPSIZE)+1" | bc)
COUNTER=0

while [ $COUNTER -lt $STEPS ]; do
  echo "Running "$KT

  #Run MC simulation at this temperature
  FILENAME=$(echo "meltingTemp_"$KT".txt")
  ./MetropolisMC $KT 5000 10000 1.0 1.0 > $FILENAME

  #Set up next iteration
  KT=$(echo "x=$KT+$STEPSIZE; if(x<1) print 0; x" | bc)
  let COUNTER=COUNTER+1
done
