#!/bin/bash
# generate vibrational constants as a function of field strength along the bond axis

for i in `seq -w -1.0 0.01 1.0`
do


k=`./6dpes BSSPV.pdb 1000.0 0.0 0.0 0.0 1000.0 0.0 0.0 0.0 1000.0 -10.0 -10.0 -10.0 0.7419 -10.0 -10.0 -10.0 0.7421 1000.0 1000.0 1000.0 1000.0 1000.0 0.0001 $i 0.0 0.0 | awk '{print $7}' | column | awk '{}{printf("%.16f\n",sqrt(($3+$1-2.0*$2)/1.0e-8)*6.811);}{}'`

echo $i $k

done

