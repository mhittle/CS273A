import tabix
import os
import sys

# REQUIRED#
# pytabix
# files containing SNP data must be in a folder specified by global variable SNP_DIR
# FILE_END global variable must be correct


# INPUT:
# this program needs 3 arguments:
#	arg1: chromosome number (1-22) or X or Y
#	arg2: start position on chromosome
#	arg3: window size to search for variants
# OUTPUT:
# run_main(): finds variants within range and returns vectors for variants
#	uses getVariants() and getVectors()
# getVariants(): returns records for each variant in range
# getVectors(): returns vectors for each record



# try:
# 'python3 getVariants.py 1 505000 50'
#	should return 6 vectors
# or:
# 'python3 getVariants.py 10 10660 50'
#	returns 32 vectors



# global variables for data file locations and names
# SNP data files with their .tbi files need to be in SNP_DIR
# all SNP files must end with FILE_END
SNP_DIR = "data/" 	# directory where SNP .gz and .gz.tbi data files are stored
FILE_END = "_ucsc_snp151.vcf.gz"	# filename ending for .gz and .gz.tbi data files




# getIdx returns the index that corresponds to the matching base
# used by getVectors() only
def getIdx(base, offset):
	bases = ['A','T','C','G']
	for i in range(len(bases)):
		if base == bases[i]:
			return offset+i
	return 0
# getFuncIDx returns the index for the corresponding vectors		
# used by getVectors() only
def getFuncIdx(func, offset):
	funcs = ['cds-indel', 'splice-3', 'splice-5', 
		'untranslated-3', 'untranslated-5', 'coding-synon', 
		'frameshift', 'intron', 'missense', 'stop-loss', 
		'ncRNA', 'nonsense', 'near-gene-3', 'near-gene-5']

	for i in range(len(funcs)):
		if func == funcs[i]:
			return offset+i
	return 0

# create vector will turn the records into 0/1 vectors
# vector id 23-dimensional with 0/1 entries
# positions are as follows:
#	1: strand of variant (0=negative, 1=positive)
#	2-5: base on reference strand (2=A, 3=T, 4=C, 5=G)
#	6-9: variant base(s) (6=A, 7=T, 8=C, 9=G)
#	10-24: variant function
#		(10:cds-indels, 11:splice-3, 12:splice-5, 13:untranslated-3,
#		14:untranslated-5, 15:coding-synonymous, 16:frameshift,
#		17:intron, 18:missense, 19:stop-loss, 20:ncRNA,
#		21:nonsense, 22:near-gene-3, 23:near-gene-5i)
# functional labels come from UCSC, more info found at https://genome.ucsc.edu/cgi-bin/hgTables
# needs getIdx() and getFuncIdx()
def getVectors(records):
	vectors = []
	n = 24	# number of fields in vector


	# turn each variant record into a vector
	for record in records:
		curr_vec = [0] * n

		# strand = plus or minus
		if record[2] == 'plus':
			curr_vec[0] = 1
		# reference base
		ref = getIdx(record[3], 2) # ref starts at idx 2
		if ref > 0:
			curr_vec[ref] = 1
		# variant bases
		variants = record[4]
		variants = variants.split('/')
		for base in variants:
			idx = getIdx(base, 6)
			curr_vec[idx] = 1
		# variant functions
		funcs = record[6]
		funcs = funcs.split(',')
		for func in funcs:
			idx = getFuncIdx(func, 10)
			curr_vec[idx] = 1



		vectors.append(curr_vec)


	## uncomment to see what vectors look like
	"""
	test = [0] * n 	# index header for visual clarity
	for i in range(n):
		test[i] = i%10
	print(test)
	for vec in vectors:
		print(vec)
	"""

	return vectors




# reads the SNP files to get all records of variants
# that fall within window_size of the start position
# RETURNS: records for variants that fall within specified window
# RECORD has 10 fields:
#	1. chromosome
#	2. start position of variant
#	3. strand of variant
#	4. refUCSC - base found on UCSC reference sequence
#	5. observed - alternate base observed
#	6. molType - molecule type, sample used to find this variant (genomic, cDNA, unknown)
#	7. Function - functional category of SNP
#	8. AllelFreqCount - number of observed alleles with frequency data
#	9. Alleles - observed alleles for which frequency data are available
#	10. allelFreqs - allele frequencies
# data comes from UCSC table browser, table snp151
# table schema can be found: https://genome.ucsc.edu/cgi-bin/hgTables
def getVariants(chrom, start_pos, window_size):
	directory = SNP_DIR 	# directory where .gz and .gz.tbi files are stored
	fn = chrom + FILE_END
	readFile = os.path.join(directory, fn)

	# make sure inputs are integers
	start_pos = int(start_pos)
	window_size = int(window_size)-1

	# open the tabix
	tb = tabix.open(readFile)

	# query for the position
	end_pos = start_pos + window_size

	print("grabbing variants from {} {} {}".format(chrom, str(start_pos), str(end_pos)))

	# grab the variant records that fall between start_pos and end_pos
	tb_records = tb.query(chrom, start_pos, end_pos)

	# store tabix data into list
	records = []
	for record in tb_records:
		records.append(record)

	return records



## use this to get an idea of how the program runs
def run_main():
	num_args = 3
	if len(sys.argv) < 1 + num_args:
		print(str(num_args) + " arguments needed.\n"
			"arg1: chromosome\n"
			"arg2: start position\n"
			"arg3: window size"
			)
		exit()


	chrom = "chr" + sys.argv[1]
	start_pos = sys.argv[2]
	size = sys.argv[3]

	# grab records for all variants in the specified window
	records = getVariants(chrom, start_pos, size) 	

	## uncomment to see what records look like
	#for record in records:
	#	print(record)

	# get vectors for all of the variant records
	vectors = getVectors(records)

	return vectors


if __name__ == "__main__":
	run_main()
