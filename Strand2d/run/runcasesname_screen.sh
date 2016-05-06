#!/bin/bash -il
# This script will run the strand2d.exec case with a StrandGen2d.input and Strand2d.input files with a user input filename saved as $1
# the output files will be saved under ./Cases/user_input_filename
screen -d -m -S Strand2dscreen
mkdir Cases/$1
cp StrandGen2d.input Strand2d.input Cases/$1
./Strand2d.exec
cp convergence.dat Cases/$1
mv StrandParts.pvtu StrandPart0.vtu output.0 solution0.pvtu residual0.pvtu surf_data.dat error.dat error0.pvtu Cases/$1
spd-say -t child_female -p 30 "Strand2d solution finished"
notify-send "Strand2d solution finished"
