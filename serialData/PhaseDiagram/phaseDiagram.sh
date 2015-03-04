#!/bin/bash

#Written by Stefan Seritan on 03/02/15 for CS 140
#Quick script that runs the MetropolisMC simulation to get phase diagram information
#PD calculated from composition of different phases as t changes

KT=0.5
STEPSIZE=0.01
MAX=1.05

STEPS=$(echo "(($MAX - $KT)/$STEPSIZE)+1" | bc)
COUNTER=0

while [ $COUNTER -lt $STEPS ]; do
  echo "Running "$KT

  #Run MC simulation at this temperature
  FILENAME=$(echo "phaseDiagram_"$KT".txt")
  ./MetropolisMC $KT 10000 50000 0.5 0.5 > $FILENAME

  #Set up next iteration
  KT=$(echo "x=$KT+$STEPSIZE; if(x<1) print 0; x" | bc)
  let COUNTER=COUNTER+1
done
