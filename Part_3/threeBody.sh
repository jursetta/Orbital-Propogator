#!/bin/bash
#
#
#

directory="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"'/'

cd $directory

clearance=(0 10 100 1000 5000 10000 50000 100000)
objective=(1 2)
tolerance=.5

for c in ${clearance[*]}; do
	for i in ${objective[*]}; do
		./threeBody $i $c $tolerance
	done;
done