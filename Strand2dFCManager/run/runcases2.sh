#!/bin/bash -il

#~ screen -d -m -S Limiter ./runcases.sh $1
#~ screen -d -m -S noLimiter ./runcases.sh $2
echo "./runcasesname.sh $1"
echo "./runcasesname.sh $2"
read -p "hit enter to start the Limiter case"
screen -d -m -S Limiter ./runcasesname.sh $1
screen -r Limiter
read -p "hit enter to start the noLimiter case"
screen -d -m -S noLimiter ./runcasesname.sh $2
screen -r noLimiter
htop
