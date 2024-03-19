#! /usr/bin/env bash


if [ "$1" == "D" ]; then
	echo "Data download"
	./launch.sh 1
elif [[ "$1" == "1" ]]; then
	echo "Performing analysis for $2 in $3... Please wait."
	./get_from_onto.sh 1 $2 $3
	./launch.sh 1a $2 $3
	./launch.sh 1b $3
	./launch.sh 2
fi
