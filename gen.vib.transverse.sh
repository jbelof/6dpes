#!/bin/bash
# generate vibrational constants as a function of field strength along the bond axis

for i in `seq -w -1.3 0.01 1.3`
do

k=`./6dpes BSSP.pdb 1000.0 0.0 0.0 0.0 1000.0 0.0 0.0 0.0 1000.0 -10.0 -10.0 -10.0 0.7419 -10.0 -10.0 -10.0 0.7421 1000.0 1000.0 1000.0 1000.0 1000.0 0.0001 0.0 0.0 $i | awk '{print $7}' | column | awk '{}{printf("%.16f\n",sqrt(($3+$1-2.0*$2)/1.0e-8)*6.811);}{}'`

echo $i $k

done

