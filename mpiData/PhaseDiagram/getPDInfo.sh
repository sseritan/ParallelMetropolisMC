#!/bin/bash

#Written by Stefan Seritan on August 26th, 2014
#Modified by Stefan Seritan on March 2nd, 2015 for CS 140 Final Project
#Get phase diagram info from MetropolisMC runs

#Control variables (initial value of kT, number of steps, step size, and counter in for loop)
KT=0.5
STEPSIZE=0.05
MAX=1.05

STEPS=$(echo "(($MAX - $KT)/$STEPSIZE)+1" | bc)
COUNTER=0

OUTFILE=pdInfo.txt

#Main for loop
while [ $COUNTER -lt $STEPS ]; do
  #Get name of job
  FILENAME=$(echo "phaseDiagram_"$KT".txt")

  #Search file for "XA (A-rich):" line
  LINE=$(grep "X1 (1-rich):" $FILENAME)
  echo $KT $LINE >> $OUTFILE

	#Get ready for next iteration
	KT=$(echo "x=$KT+$STEPSIZE; if(x<1) print 0; x" | bc)
	let COUNTER=COUNTER+1
done
