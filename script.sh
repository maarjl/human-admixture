#!/bin/sh

## HUMAN ADMIXTURE PIPELINE : Beagle v4 -> ChromoPainter v2 -> GLOBETROTTER
## Author: Maarja Lepamets
## 
## Institute of Molecular and Cellular Biology
## University of Tartu, Estonia
##
## January, 2015


## constants
DATAPATH="./data/"
PREFIX="xingnorel"

echo
echo "Starting HUMAN ADMIXTURE PIPELINE"
echo "Searching data from path: $DATAPATH"
echo "Prefix of data files: $PREFIX"
echo

## .bed to .ped
if [ ! -f ${DATAPATH}${PREFIX}.ped ]; then 
	printf "Converting ${DATAPATH}${PREFIX}.bed to PED format..."
	time -f %E ./plink --noweb --bfile ${DATAPATH}${PREFIX} --recode --tab --out ${DATAPATH}${PREFIX} --silent
	echo
	# Takes 1.5 min
fi

## phasing pipeline
if [ -f ${DATAPATH}${PREFIX}.phased.vcf.gz ]; then
	gunzip ${DATAPATH}${PREFIX}.phased.vcf.gz
fi
if [ ! -f ${DATAPATH}${PREFIX}.phased.vcf ]; then
	echo "Phasing the data..."
	mkdir -p ${DATAPATH}chrom
	# change the loop if you have different chromosome numbers than 1 to 22
	for i in $(seq 1 22); do
		./pipeline/beagle.sh ${DATAPATH} ${PREFIX} $i &
	done
	wait
	echo	
	printf "Concatenating data..."
	time -f %E ./pipeline/union.sh ${DATAPATH}chrom/${PREFIX} > ${DATAPATH}${PREFIX}.phased.vcf
	echo
	gzip ${DATAPATH}${PREFIX}.phased.vcf
	#rm -fr ${DATAPATH}chrom
fi

## Beagle output to ChromoPainter input conversion
if [ ! -f ${DATAPATH}${PREFIX}.haplotypes ] || [  ! -f ${DATAPATH}${PREFIX}.recomrates ]; then
	printf "Creating input files for ChromoPainter v2..."
	time -f %E python ./pipeline/beagle_to_chromopainter_convert.py ${DATAPATH}${PREFIX}
fi
nsamples=$(./pipeline/create_population_list_infile_and_idfile.py ${DATAPATH}${PREFIX} $@)

## ChromoPainter
printf "Running ChromoPainterv2..."
time -f %E ./pipeline/chromopainter.sh ${DATAPATH} ${PREFIX} $nsamples

## GLOBETROTTER
printf "Running GLOBETROTTER..."
time -f %E R < GLOBETROTTER.R ${DATAPATH}${PREFIX}.param ${DATAPATH}${PREFIX}.samples.out ${DATAPATH}${PREFIX}.recomrate --no-save > ${DATAPATH}${PREFIX}.globetrotter.log










