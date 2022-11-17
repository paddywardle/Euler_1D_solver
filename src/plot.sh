#!/bin/bash

file="../data/test_LF.dat"

for i in {2..5}
do
	echo "plot '$file' using 1:$i with linespoints" | gnuplot -persist
done	
