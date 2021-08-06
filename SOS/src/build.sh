#!/bin/bash

program=MC
gfortran-mp-4.8 Public_MC.f90 MonteCarloRayleigh_SOS.f90 -o $program
rm -rf *.mod

num=`cat Count`
echo $num
dir=./output_$num
mkdir $dir

 ./$program
mv trajectory.dat radiance.dat  $dir

echo "data stored in directory" $num
num=$(($num+1))
echo $num > ./Count
