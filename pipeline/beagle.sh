#!/bin/sh

## Beagle pipeline for one chromosome

# current chromosome number
chromosome=$3

# get start of current chromosome and number of SNPs from current chromosome
first=$(($(grep -P -n "^${chromosome}\t" ${1}${2}.map | cut -f 1 -d : | head -n 1)+6))
last=$(($(grep -P -n "^${chromosome}\t" ${1}${2}.map | cut -f 1 -d : | tail -n 1)+6))


# .ped to .vcf

# temporary files
cut -f 1-6,${first}-${last} ${1}${2}.ped > ${1}chrom/${2}.$3.ped
grep -P "^${chromosome}\t" ${1}${2}.map > ${1}chrom/${2}.$3.map

# convert to .vcf
./pipeline/bioinformatics_format_convert.py input_file_1=${1}chrom/${2}.$3.ped input_file_2=${1}chrom/${2}.$3.map input_type=PLINK output_file_1=${1}chrom/${2}.$3.vcf output_type=VCF
# phase
java -Xmx4000m -jar ./beagle.r1398.jar gt=${1}chrom/${2}.$3.vcf out=${1}chrom/${2}.$3.phased nthreads=4 > ${1}chrom/log$3
echo "Finished phasing chromosome $3..."

