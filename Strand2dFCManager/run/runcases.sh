#!/bin/bash -il
# This script will run the Strand2dFC case with a user input input.namelist
#~ read -p "Please input your namelist " N
read -p "please input your namelist  " N
mkdir Cases/$N
cp $N Cases/$N
./Strand2dFCSolver.exec $N
mv output.0 strandSolution0.pvtu error.dat forces.dat pressure.dat Cases/$N
