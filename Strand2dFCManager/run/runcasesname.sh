#!/bin/bash -il
# This script will run the Strand2dFC case with a user input input.namelist
mkdir Cases/$1
cp $1 Cases/$1
./Strand2dFCSolver.exec $1
mv output.0 strandSolution0.pvtu error.dat forces.dat pressure.dat Cases/$1
cp convergence.dat Cases/$1
spd-say -t child_female -p 30 "Strand2dFC solution finished"
notify-send "Strand2dFC solution finished"
