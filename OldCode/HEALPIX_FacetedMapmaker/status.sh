#!/bin/sh

echo -e "\n--------------------------"

for ((  comp = 0 ;  comp <= 14;  comp++  ))
do
	if [ $comp -ge 10 ]
	then
		echo -n "On eor-$comp: "
		echo -n "Processes running: "
		ssh eor-$comp "ps -C facetedMapmaker | grep -c facet"
	else
		echo -n "On eor-0$comp: "
		echo -n "Processes running: "
		ssh eor-0$comp "ps -C facetedMapmaker | grep -c facet"
	fi
done
