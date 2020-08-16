#!/bin/bash

#rootfiles="/data/11b/itabashi/root/*"
rootfiles="/data/41a/ELS/okuyama/rootfiles/nnL_small/*"

for pathfile in $rootfiles; do
	echo $pathfile >> "small2.list"
done
	
