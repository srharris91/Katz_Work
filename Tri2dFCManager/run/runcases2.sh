#!/bin/bash -il
# this will run 2 simultaneous cases of runcasename.sh in 2 screen values using the 2 input.namelist values 
# to run use the cmd line "./runcases2.sh input1.namelist and input2.namelist" where the input.namelist's are the names of your specific namelists
echo "./runcasesname.sh $1"
echo "./runcasesname.sh $2"
read -p "hit enter to start the $1 case"
screen -d -m -S $1 ./runcasesname.sh $1
screen -r $1
read -p "hit enter to start the $2 case"
screen -d -m -S $2 ./runcasesname.sh $2
screen -r $2
htop
