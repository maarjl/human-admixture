#!/bin/bash

mkdir -p $1EMest
# number of samples for parameter estimates
div=10
s=$(($3 / $div))
if [ $s -le 5 ]; then
	s=5
fi

# estimating parameters
for i in $(seq 1 $s); do
	r=$RANDOM
	n=$(( $r % $3 ))
	../ChromoPainterv2 -a 0 0 -i 10 -in -iM -s 0 -g ${1}${2}.haplotypes -r ${1}${2}.recomrates -t ${1}${2}.idfile -f ${1}${2}.poplist $n $n -o $1EMest/${2}.$n > $1EMest/log.$n &
done
wait

# calculating final parameters for the actual run
./neaverage.pl -o ${1}${2}.neaverage.txt ${1}EMest/${2}.*.EMprobs.out
necmd=$(cat ${1}${2}.neaverage.txt)

# actual ChromoPainter run
../ChromoPainterv2 -s 10 $necmd -g ${1}${2}.haplotypes -r ${1}${2}.recomrates -t ${1}${2}.idfile -f ${1}${2}.poplist 0 0 -o ${1}${2} > ${1}${2}.log









