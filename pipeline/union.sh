#!/bin/sh

# Makes union of all chromosome files

path=$1

for i in $(seq 21 22); do
	gunzip $1.$i.phased.vcf.gz
	if [ $i -eq 21 ]; 
	then
		# get the whole first chromosome file
		cat $1.$i.phased.vcf
	else
		# get data from chromosome file
		grep -v "#" $1.$i.phased.vcf
	fi
done