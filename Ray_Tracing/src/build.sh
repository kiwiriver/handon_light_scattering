#!/bin/bash
program=GOM2D

gfortran-mp-4.8 Public_GOM.f90 Public_Sphere.f90 GOM2D.f90 -o $program
rm -rf *.mod

num=`cat Count`
echo $num
dir=./output_$num 
mkdir $dir

./$program
mv log.dat path.dat p11.dat shape.dat  $dir

echo "data stored in directory" $num
num=$(($num+1))
echo $num > ./Count
