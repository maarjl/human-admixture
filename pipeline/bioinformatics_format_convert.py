#!/usr/bin/python

"""
Copyright (c) 2012, 2013 The PyPedia Project, http://www.pypedia.com
<br>All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met: 

# Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer. 
# Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution. 

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

http://www.opensource.org/licenses/BSD-2-Clause
"""

__pypdoc__ = """
Method: bioinformatics_format_convert
Link: http://www.pypedia.com/index.php/bioinformatics_format_convert
Retrieve date: Sat, 24 Jan 2015 12:58:19 +0200



Convert between various bioinformatics formats. The supported formats are:
* [http://pngu.mgh.harvard.edu/~purcell/plink/tutorial.shtml PLINK] ped and map format
* [http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#tr TPLINK] tped and tfam format
* [http://faculty.washington.edu/browning/beagle/beagle_3.3.2_31Oct11.pdf  BEAGLE] bgl and markers format
* [http://www.stats.ox.ac.uk/~marchini/software/gwas/file_format.html IMPUTE2] gen and sample format
* [http://www.sph.umich.edu/csg/abecasis/Merlin/tour/input_files.html MERLIN] ped and dat format
* [http://www.1000genomes.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-41 VCF] format

'''Options''':
* input_file_1 : The first input file
** This should be the PED file for PLINK, the TPED file for TPLINK, the beagle file for BEAGLE, the GEN file for IMPUTE2, the PED file for MERLIN and the VCF file for VCF.
* input_file_2 : The second input file
** This should be the MAP file for PLINK, the TFAM file for TPLINK, the markers file for BEAGLE, the sample file IMPUTE2, the DAT file for MERLIN and None for VCF.
* input_type : Available values are: PLINK, TPLINK, BEAGLE, IMPUTE2, MERLIN, VCF
* output_file_1, output_file_2 the correspondent output files according to selected output format. For VCF this can be None.
* output_type : The format of the output files. The available options are the same as with input_type
* chromosome : In case the input format does not have chromosome information (i.e. BEAGLE) you can define it here. Only string is allowed
* phenotype : In case the input format supports multiple phenotypes (i.e. IMPUTE2) you can define which one should be picked.
* gender : In case the input format supported multiple phenotypes and does not make any distinction between a regular phenotype and gender (i.e. IMPUTE2) you can put the name of the phenotype that corresponds to gender.
* silent: set True to suppress output

[[Category:Validated]]
[[Category:Algorithms]]
[[Category:Bioinformatics]]
[[Category:Format Conversion]]


"""

import os
import re
import glob
import gzip
import numpy
import tempfile
import mimetypes
import itertools
import sys


class progress_bar():
	def __init__(self, iterations):
		self.iterations = iterations
		self.prog_bar = '[]'
		self.fill_char = '*'
		self.width = 40
		self.__update_amount(0)
		self.animate = self.animate_noipython

	def update_iteration(self, elapsed_iter):
		self.__update_amount((elapsed_iter / float(self.iterations)) * 100.0)
		self.prog_bar += '  %d of %s complete' % (elapsed_iter, self.iterations)

	def __update_amount(self, new_amount):
		percent_done = int(round((new_amount / 100.0) * 100.0))
		all_full = self.width - 2
		num_hashes = int(round((percent_done / 100.0) * all_full))
		self.prog_bar = '[' + self.fill_char * num_hashes + ' ' * (all_full - num_hashes) + ']'
		pct_place = (len(self.prog_bar) // 2) - len(str(percent_done))
		pct_string = '%d%%' % percent_done
		self.prog_bar = self.prog_bar[0:pct_place] + \
			(pct_string + self.prog_bar[pct_place + len(pct_string):])

	def __str__(self):
		return str(self.prog_bar)

class bioinformatics_file_helper:
	'''
	This is a collection of common functions used commonly
	when opening and processing files with genotype information
	'''

	@staticmethod
	def get_alleles(genotypes):
		'''
		genotypes: a list with tuples of genotypes for the same SNP. i.e: [('A', 'A'), ('A', 'G'), ('G', '0')]
		returns: a tuple with the different alleles. i.e ('A', 'G')
		The alleles in the tuple are sorted according to their frequency.
		'''

		flat_genotypes = [genotype for genotype_pair in genotypes for genotype in genotype_pair]
		all_alleles = list(set(flat_genotypes) - set(['0']))
		all_alleles_c = len(all_alleles)

		if all_alleles_c == 0:
			return ('0', '0')
		if all_alleles_c == 1:
			return (all_alleles[0], '0')
		else:
			count_sorted = {x : flat_genotypes.count(x) for x in all_alleles}
			return sorted(count_sorted, key=lambda x : count_sorted[x])[::-1]  # Sort and revert

	@staticmethod
	def get_chromosome_files(path, chromosome_exp=r'chr%(chromosome)s'):
		'''
		Check if there are files that match 'chromosome_exp'
		Check for multiple matches. 
		This method is useful for getting genomic files split per chromosome

		Return: A pattern that match the files and the chromosomes that matched. For example if the files are:
		path/data_chr1.vcf
		path/data_chr2.vcf

		The return value of get_chromosome_files('path/*.vcf') will be:

		('path/data_chr%(chromosome)s.vcf', [1,2])
		'''

		reference_files = glob.glob(path)
		if not reference_files:
			# print 'Warning: path %s is empty' % path
			return None, None
		reference_filenames = [os.path.split(x)[1] for x in reference_files]

		if '%(chromosome)s' not in chromosome_exp:
			print "Warning: Could not find chromosome identifier '%(chromosome)s' in parameter chromosome_exp"
			return None, None

		chr_reg_expr = chromosome_exp.replace('%(chromosome)s', '([\d]+|X|Y)')

		# Check if there is a chromosome identifier in the filenames
		stem_d = {}
		chromosomes = []
		for reference_filename in reference_filenames:
			s = re.search(chr_reg_expr, reference_filename)
			if s:
				chromosome = s.group(1)
				if chromosome in [str(x) for x in range(1, 23)] + ['X', 'Y']:
					stem = reference_filename.replace(s.group(0), chromosome_exp)
					if stem_d.has_key(chromosome):
						print 'Warning: Do not know which of these two files is the proper file for chromosome %s: %s, %s' % (chromosome, stem_d[chromosome] % {'chromosome' : chromosome}, reference_filename)
						return None, None
					else:
						stem_d[chromosome] = stem
						chromosomes += [chromosome]
				else:
					print 'Warning: Ignoring file: %s . Cannot recognize chromosome: %s' % (reference_filename, chromosome)
		# Check if all stems are the same:
		stem_set = set([values for values in stem_d.itervalues()])
		if len(stem_set) > 1:
			print 'Warning: the \'%s\' part should be in the same position for all the files' % (chromosome_exp)
			print 'These are not the same:'
			print '\n'.join(stem_set)
			return None, None
		if len(stem_set) == 0:
			print 'Warning: could not find %s in any file in path %s' % (chromosome_exp, path)
			return None, None

		return list(stem_set)[0], chromosomes

	@staticmethod
	def line_reader(f):
		'''
		f: an opened file 
		Returns a list of identifiers in a line that is whitespace separated
		'''

		return f.readline().replace('\n', '').split()

	@staticmethod
	def line_writer(f, line, sep='\t'):
		'''
		f: a file
		line: a list of values

		Saves the list 'line' in file 'f' separated with 'sep'
		'''
		return f.write(sep.join(line) + '\n')

	@staticmethod
	def open_file_read(filename, mode='rU'):
		'''
		Checks the type of filename and returns a file or a string stream
		'''
		if type(filename) is str:
			# Check if file is a gzip
			if mimetypes.guess_type(filename)[1] == 'gzip':
				return gzip.open(filename, 'rb')

			return open(filename, mode)
		else:
			filename.seek(0)
			return filename

	@staticmethod
	def open_file_write(filename):
		'''
		Checks the type of filename and returns a file or a string stream
		'''
		if type(filename) is str:

			if mimetypes.guess_type(filename)[1] == 'gzip':
				return gzip.open(filename, 'wb')

			return open(filename, 'w')
		else:
			return filename

	@staticmethod
	def close_file(stream):
		'''
		Checks the type of filesource and closes the file
		'''
		if type(stream) is file:
			stream.close()

	@staticmethod
	def line_generator(filename):
		'''
		filename: a filename or open file
		Generates lines from a filename. Lines are splitted on whitespaces
		'''

		read_from = bioinformatics_file_helper.open_file_read(filename)

		line_counter = 0
		for l in read_from:
			line_counter += 1
			s = l.replace('\n', '').split()
			yield line_counter, s

		bioinformatics_file_helper.close_file(read_from)

	@staticmethod
	def column_generator(filename, batch_size=10000):
		'''
		filename: a filename or open file
		Reads a column of a file.
		After 'batch_size' reads reopens the file
		yields a tuple: current column, line
		'''

		start_column = 0
		line_counter = 0
		while True:
			to_return = []
			finished = False

			f = bioinformatics_file_helper.open_file_read(filename)
			for l in f:
				s = l.replace('\n', '').split()
				# Store a batch_size at most records
				current_line_batch = s[start_column: start_column + batch_size]
				
				if not len(current_line_batch):
					finished = True
					break

				to_return += [current_line_batch]

			if finished:
				break

			bioinformatics_file_helper.close_file(f)

			# Transpose
			to_return_transposed = numpy.transpose(to_return)

			# Yield columns
			for line in to_return_transposed:
				line_counter += 1
				yield line_counter, list(line)

			start_column += batch_size

	@staticmethod
	def column_writer(filename, batch_size=10000, silent=False):
		'''
		Saves to file column by column.
		It is implemented as a generator
		If this method consumes too much memory, try lowering the batch_size    
		To suppress output set silent=True

		Example:
		g = write_file_vertical('test.txt')
		g.next()
		g.send(['1', '2'])
		g.send(['3', '4'])
		g.send(['5', '6'])
		g.send(None)
		'''

		finished = False
		while not finished:
			current_batch = []
			for current_record in range(batch_size):
				data = (yield True)
				if not data:
					finished = True
					break
				current_batch += [data]

			if current_batch:
				current_batch_transposed = numpy.transpose(current_batch)

				new_temp_file = tempfile.NamedTemporaryFile(delete=False)
				if not silent:
					print 'Created: ', new_temp_file.name

				for line in current_batch_transposed:
					bioinformatics_file_helper.line_writer(new_temp_file, line)

				old_temp_filename = new_temp_file.name
				new_temp_file.close()

		if type(filename) is str:
			os.rename(old_temp_filename, filename)
			if not silent:
				print 'Moved %s to %s' % (old_temp_filename, filename)
		else:
			# Dump temporary file to stream. Then delete it
			temp_file = bioinformatics_file_helper.open_file_read(old_temp_filename)
			for l in temp_file:
				filename.write(l)
			bioinformatics_file_helper.close_file(old_temp_filename)
			os.unlink(old_temp_filename)
			if not silent:
				print 'Delete: ', old_temp_filename
		yield False




def BEAGLE_writer(beagle_filename, markers_filename, reader):
	'''
	Write in beagle genotype format
	Description: http://faculty.washington.edu/browning/beagle/beagle_3.3.2_31Oct11.pdf 
	'''

	bfh = bioinformatics_file_helper()

	beagle_file = bfh.open_file_write(beagle_filename)
	markers_file = bfh.open_file_write(markers_filename)

	header = reader.next()
	sample_names = header['sample_ids']
	phenotype_ids = header['phenotype_ids']

	identifier_line = ['I', 'id'] + [sample_name for sample_name_pair in zip(sample_names, sample_names) for sample_name in sample_name_pair]
	bfh.line_writer(beagle_file, identifier_line)

	affection_line = [
		'A',
		header['phenotype_name'] if header.has_key('phenotype_name') else 'phenotype',
		] + [phenotype_id for phenotype_id_pair in zip(phenotype_ids, phenotype_ids) for phenotype_id in phenotype_id_pair]
	bfh.line_writer(beagle_file, affection_line)

	for record in reader:

		genotypes = record['genotypes']

		markers_line = [
			record['rs_id'],
			record['position'],
		] + list(bfh.get_alleles(genotypes))
		bfh.line_writer(markers_file, markers_line)

		beagle_line = [
			'M',
			record['rs_id'],
		] + [g for g_pair in genotypes for g in g_pair]
		bfh.line_writer(beagle_file, beagle_line)

	bfh.close_file(beagle_file)
	bfh.close_file(markers_file)


def PLINK_writer(ped_filename, map_filename, reader, genotypes_per_batch=10000, silent=False):
	'''
	Save plink's format ped_map file
	format description: http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#ped
	'''

	bfh = bioinformatics_file_helper()

	# Get header
	header = reader.next()

	# Save header
	ped_writer = bfh.column_writer(ped_filename, genotypes_per_batch, silent=silent)
	ped_writer.next()
	ped_writer.send(header['family_ids'])
	ped_writer.send(header['sample_ids'])
	ped_writer.send(header['father_ids'])
	ped_writer.send(header['mother_ids'])
	ped_writer.send(header['sex_ids'])
	ped_writer.send(header['phenotype_ids'])

	# Save markers
	map_file = bfh.open_file_write(map_filename)
	for marker in reader:
		ped_writer.send([genotype[0] for genotype in marker['genotypes']])
		ped_writer.send([genotype[1] for genotype in marker['genotypes']])

		map_record = [
			marker['chromosome'],
			marker['rs_id'],
			'0',
			marker['position']
		]
		bfh.line_writer(map_file, map_record)

	# Close ped file
	ped_writer.send(None)

	# Close map file
	bfh.close_file(map_file)

def VCF_writer(vcf_filename, reader):
	'''
	Description: http://www.1000genomes.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-41 
	'''

	bfh = bioinformatics_file_helper()

	def vcf_genotype(ref, alt, genotype):
		if genotype == '0':
			return '.'
		elif genotype == ref:
			return '0'
		elif genotype == alt:
			return '1'
		else:
			raise Exception('more than one alternatives are not supported by this converter')

	vcf_file = bfh.open_file_write(vcf_filename)

	header = reader.next()

	vcf_header = [['##fileformat=VCFv4.1']]
	vcf_header += [['##source=pypedia.com']]
	vcf_header += [['##FILTER=<ID=PASS,Description="Passed variant FILTERs">']]
	vcf_header += [['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT'] + header['sample_ids']]

	#Write header
	for entry in vcf_header:
		bfh.line_writer(vcf_file, entry)

	for record in reader:
	
		alleles = bfh.get_alleles(record['genotypes'])

		to_write = [
			record['chromosome'],
			record['position'],
			record['rs_id'],
			alleles[0],
			alleles[1] if alleles[1] <> '0' else 'N',
			'.', '.', '.',
			'GT',
		] +   ['/'.join([vcf_genotype(alleles[0], alleles[1], allele) for allele in gen_pair]) for gen_pair in record['genotypes']]
		
		bfh.line_writer(vcf_file, to_write)

	bfh.close_file(vcf_file)
	

def VCF_reader(vcf_filename):
	'''
	Description: http://www.1000genomes.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-41
	'''

	bfh = bioinformatics_file_helper()

	def vcf_split(vcf_gen):
		if '/' in vcf_gen:
			tmp_split = vcf_gen.split('/') 
		elif '|' in vcf_gen:
			tmp_split = vcf_gen.split('|')
		else:
			raise Exception('Invalid field: %s in genotypes' % vcf_gen)

		return [0 if x == '.' else (int(x)+1) for x in tmp_split]

	vcf_generator = bfh.line_generator(vcf_filename)

	for vcf_c, vcf_s in vcf_generator:
		if vcf_s[0][0:2] == '##':
			continue

		if vcf_s[0] == '#CHROM':
			sample_ids = vcf_s[9:]
			samples = len(sample_ids)
			yield {
				'sample_ids' : sample_ids,
				'family_ids' : ['1'] * samples,
				'father_ids' : ['0'] * samples,
				'mother_ids' : ['0'] * samples,
				'sex_ids' : ['0'] * samples,
				'phenotype_ids' : ['-9'] * samples,
			}
			continue

		ref_allele = vcf_s[3]
		alt_allele = vcf_s[4].split(',')
		all_alleles = ['0', ref_allele] + alt_allele
		format = vcf_s[8]
		genotype_index = format.split(':').index('GT')

		chromosome = vcf_s[0]
		if 'chr' in chromosome.lower():
			chromosome = chromosome[3:]

		yield {
			'chromosome' : chromosome,
			'position' : vcf_s[1],
			'rs_id' : vcf_s[2],
			'genotypes' : [tuple([all_alleles[ar_allele] for ar_allele in vcf_split(vcf_gen.split(':')[genotype_index])]) for vcf_gen in vcf_s[9:]]
		}	


def BEAGLE_reader(beagle_filename, marker_filename, chromosome):
	'''
	Reader for beagle data
	description: http://faculty.washington.edu/browning/beagle/beagle_3.3.2_31Oct11.pdf 
	'''

	if chromosome is None:
		raise Exception('To read from beagle files, "chromosome" parameter should be set')

	bfh = bioinformatics_file_helper()

	# Get header 
	identifier_found, affection_found = False, False 
	for beagle_c, beagle_s in bfh.line_generator(beagle_filename):

		if beagle_s[0] == 'I':
			identifier_line = beagle_s
			identifier_found = True

		elif beagle_s[0] == 'A':
			affection_line = beagle_s
			affection_found = True

		if beagle_c > 10 and identifier_found:
			break

	if not identifier_found:
		raise Exception('Could not find Identifier (I) line in ' + beagle_filename)

	# Yield header
	samples = len(identifier_line[2::2])
	header = {
		'family_ids' : ['1' for sample in range(samples)],
		'sample_ids' : identifier_line[2::2],
		'father_ids' : ['0' for sample in range(samples)],
		'mother_ids' : ['0' for sample in range(samples)],
		'sex_ids' : ['0' for sample in range(samples)],
		'phenotype_ids' : ['0' for sample in range(samples)],
	}
	
	if affection_found:
		header['phenotype_ids'] = affection_line[2::2]
		header['phenotype_name'] = affection_line[1]

	yield header

	# Yield rest of the data
	beagle_generator = bfh.line_generator(beagle_filename)
	marker_generator = bfh.line_generator(marker_filename)

	for beagle_data, marker_data in itertools.izip(beagle_generator, marker_generator):

		# Skip header
		while beagle_data[1][0] in ['A', 'I']:
			beagle_data = beagle_generator.next()

		yield {
			'chromosome' : chromosome,
			'rs_id' : marker_data[1][0],
			'position' : marker_data[1][1],
			'genotypes' : [(beagle_data[1][i], beagle_data[1][i + 1]) for i in range(2, (samples * 2) + 1, 2)],
		}

def PLINK_reader(ped_filename, map_filename, genotypes_per_batch=10000):
	'''
	generator for a plink's PED and MAP file
	format description: http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#ped
	'''

	bfh = bioinformatics_file_helper()

	chromosomes = []
	rs_ids = []
	positions = []

	# Load mapfile
	for map_c, map_s in bfh.line_generator(map_filename):
		chromosomes += [map_s[0]]
		rs_ids += [map_s[1]]
		positions += [map_s[3]]

	# Read ped file column by column
	ped_reader = bfh.column_generator(ped_filename, genotypes_per_batch)

	# yield header 
	yield {
		'family_ids' : ped_reader.next()[1],
		'sample_ids' : ped_reader.next()[1],
		'father_ids' : ped_reader.next()[1],
		'mother_ids' : ped_reader.next()[1],
		'sex_ids' : ped_reader.next()[1],
		'phenotype_ids' : ped_reader.next()[1],
	}

	# yield genotypes
	for marker in range(len(rs_ids)):
		yield {
			'chromosome' : chromosomes[marker],
			'position' : positions[marker],
			'rs_id' : rs_ids[marker],
			'genotypes' : zip(ped_reader.next()[1], ped_reader.next()[1])
		}


def bioinformatics_format_convert(
	input_file_1,
	input_file_2,
	input_type,
	output_file_1,
	output_file_2,
	output_type,
	chromosome=None,
	phenotype='pheno',
	gender='gender',
	genotypes_per_batch=10000,
	silent=True):

	if input_type == 'PLINK':
		reader = PLINK_reader(input_file_1, input_file_2)
	elif input_type == 'BEAGLE':
		reader = BEAGLE_reader(input_file_1, input_file_2, chromosome)
	elif input_type == 'VCF':
		reader = VCF_reader(input_file_1)
	else:
		raise Exception('Unknowm file type: %s in parameter input_type' % (str(input_type)))

	if output_type == 'PLINK':
		PLINK_writer(output_file_1, output_file_2, reader, silent=silent)
	elif output_type == 'BEAGLE':
		BEAGLE_writer(output_file_1, output_file_2, reader)
	elif output_type == 'VCF':
		VCF_writer(output_file_1, reader)
	else:
		raise Exception('Unknown file type: %s in parameter output_type' % (str(output_type)))

	if not silent:
		print 'Written output file:', str(output_file_1)
		if output_file_2:
			print 'Writen output file:', str(output_file_2)



arguments = {"input_file_1":"", "input_file_2":"", "input_type":"", "output_file_1":"", "output_file_2":"", "output_type":""}
# Method name =bioinformatics_format_convert()
if __name__ == '__main__':
	for i in range(1, len(sys.argv)):
		key, value = sys.argv[i].split("=", 1)
		arguments[key] = value
	#print arguments
	#chromosome = None
	#phenotype = 'pheno'
	#gender = 'gender'
	
	returned = bioinformatics_format_convert(input_file_1 = arguments["input_file_1"], 
											input_file_2 = arguments["input_file_2"], 
											input_type = arguments["input_type"], 
											output_file_1 = arguments["output_file_1"], 
											output_file_2 = arguments["output_file_2"], 
											output_type = arguments["output_type"])
	if returned:
		print 'Method returned:'
		print str(returned)
	
		
