#!/bin/bash
for ((i=0;i<=300;i++))
do
python calc_slme.py OUTCAR 0.81 0.81  ${i}E-6  -0.49
done 
