#!/bin/sh
# can be used to change file names quickly
ls *SPLam* | while read f
do
    mv $f ${f/SPLam/SPTurbSA}
done