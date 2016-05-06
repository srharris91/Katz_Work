#!/bin/bash -il
# This script will run the Tri2dFCSolver.exec case with a user input input.namelist
mkdir Cases/$1
cp $1 Cases/$1
./Tri2dFCSolver.exec $1
cp convergence.dat Cases/$1
mv output.0 triSolution0.pvtu error.dat forces.dat Cases/$1
spd-say -t child_female -p 30 "Tri2dFC solution finished"
notify-send "Tri2dFC solution finished"
