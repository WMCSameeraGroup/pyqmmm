#!/bin/csh

set INP=waterbox.com
set SCR=wat

setenv GAUSS_SCRDIR /path/to/scratch/folder/${SCR}

mkdir -p $GAUSS_SCRDIR

g16 < ${INP} > ${INP:r}.log

rm -rf $GAUSS_SCRDIR

