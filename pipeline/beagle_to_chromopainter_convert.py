#!/usr/bin/python 

import sys

# takes info from vcb and creates haplotypes and recombination rate files
def vcf_to_haplotypes_and_recomrates_convert (filename_prefix):
    input_vcf_filename = filename_prefix + ".phased.vcf"
    output_haplotypes_filename = filename_prefix + ".haplotypes"
    output_recomrates_filename = filename_prefix + ".recomrates"
    
    input_vcf_file = open(input_vcf_filename, "r")
    
    offset = 0
    is_first_data_row = False
    for line in input_vcf_file:
        offset += len(line)
        if is_first_data_row:
            break
        if line[0:6] == "#CHROM":
            n_individuals = len(line.split("\t")) - 9
            is_first_data_row = True
    data = line.strip().split("\t")
    last_SNP = data[0]
    positions = [data[1]]
    all_SNP_data = [data[9:]]
    
    # writes recombination rate file
    output_recomrates_file = open(output_recomrates_filename, "wb")
    output_recomrates_file.write("start.pos\trecom.rate.perbp\n")
    for line in input_vcf_file:
        data = line.strip().split("\t")
        output_recomrates_file.write(positions[-1] + "\t" + ("0.0000001" if last_SNP == data[0] else "-9") + "\n")
        last_SNP = data[0]
        positions.append(data[1])
        all_SNP_data.append(data[9:])
    output_recomrates_file.write(positions[-1] + "\t0\n") 
    output_recomrates_file.close()
    
    # writes haplotypes file
    output_haplotypes_file = open(output_haplotypes_filename, "wb")
    output_haplotypes_file.write(str(n_individuals * 2) + "\n")
    output_haplotypes_file.write(str(len(positions)) + "\n")
    output_haplotypes_file.write("P " + " ".join(positions) + "\n")
    for i in range(n_individuals):
        ref = alt = ""
        for SNP_data in all_SNP_data:
            ref += SNP_data[i][0]
            alt += SNP_data[i][2]
        output_haplotypes_file.write(ref + "\n" + alt + "\n")
    output_haplotypes_file.close()
    
    input_vcf_file.close()
    

prefix = sys.argv[1]
vcf_to_haplotypes_and_recomrates_convert (prefix)