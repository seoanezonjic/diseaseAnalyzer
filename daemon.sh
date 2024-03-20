#! /usr/bin/env bash


if [ "$1" == "D" ]; then
	echo "Data download"
	./launch.sh config/test 1
elif [[ "$1" == "1" ]]; then
	echo "Performing analysis for $2 in $3... Please wait."
	./launch.sh config/test HPO $2 $3
	./launch.sh config/test 1a $2 $3
	./launch.sh config/test 1b $3
	./launch.sh config/test 2
fi
