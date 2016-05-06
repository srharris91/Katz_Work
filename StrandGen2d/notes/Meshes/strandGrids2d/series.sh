#!/bin/sh
################################################################################
# Script to determine solution error of a scheme
################################################################################


#---------------------------begin configure section----------------------------#

# number of surface cells
NPS=(4 8 16 32 64 128 256 512)

# perturbabion
PERTURB=1

# skew angle
SKEW=5

#---------------------------end configure section------------------------------#

make
N=`expr 0`
for NP in ${NPS[*]}
do
    echo "running grid" line$NP'p'$PERTURB'a'$SKEW'.strand2d'
    ./line -n ${NPS[$N]} -p $PERTURB -a $SKEW -m line$NP'p'$PERTURB'a'$SKEW'.strand2d'
    N=`expr $N + 1`
done