#!/bin/bash
# generate vibrational frequency along the minimum

xvec=-0.5572
yvec=-0.5497
zvec=-0.6223

xmin=-4.176
ymin=-4.216
zmin=8.437


for i in `seq -w -1.0 0.001 1.0`
do

x=`echo $xmin | awk -v inc=$i -v xvec=$xvec '{}{printf("%.3f", xvec*inc + $1)}{}'`
y=`echo $ymin | awk -v inc=$i -v yvec=$yvec '{}{printf("%.3f", yvec*inc + $1)}{}'`
z=`echo $zmin | awk -v inc=$i -v zvec=$zvec '{}{printf("%.3f", zvec*inc + $1)}{}'`

k=`./6dpes MOF5+H2.globalmin.pdb 25.669 0.0 0.0 0.0 25.669 0.0 0.0 0.0 25.669 $x $y $z 0.7419 $x $y $z 0.7421 1000.0 1000.0 1000.0 1000.0 1000.0 0.0001 0.0 0.0 0.0 | awk '{print $7}' | column | awk '{}{printf("%.16f\n",sqrt(($3+$1-2.0*$2)/1.0e-8)*42.796);}{}'`

echo $i $k


done


