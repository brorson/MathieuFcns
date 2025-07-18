#!/bin/bash

rm mathieu*.csv
rm ../MathieuStudies/mathieu*.csv

fns=(-1000 -100 -10 -1 -0.1 -0.001 0.001 0.01 0.1 1 10 100 1000)

for i in ${fns[@]};
do
echo Writing ce values to file for q = $i ...
echo $i | octave write_mathieu_ce_gvs.m
done
echo ================================================

for i in ${fns[@]};
do
echo Writing se values to file for q = $i ...    
echo $i | octave write_mathieu_se_gvs.m
done
echo ================================================


for i in ${fns[@]};
do
echo Writing ce deriv values to file for q = $i ...
echo $i | octave write_mathieu_ce_deriv_gvs.m
done
echo ================================================


for i in ${fns[@]};
do
echo Writing se deriv values to file for q = $i ...
echo $i | octave write_mathieu_se_deriv_gvs.m
done
echo ================================================


echo Copying files over
cp *.csv ../MathieuStudies/

