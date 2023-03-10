#!/usr/bin/env python

'''Summarize average genome size and genome equivalents from output.

This script takes the following input:
- directory containing *.census files

This script produces the following output:
- tsv file with microbe census data from each sample

'''

import argparse, csv, glob, os
from collections import defaultdict
import pandas as pd


def parse_census(input_directory, verbose):
	# generate dictionary of gene set coverage & breadth values
	# initialize dictionary & list objects

	os.chdir(input_directory)					# change to input directory
	file_list = glob.glob("*.census")			# generate list of *.census files

	file_counter = 0
	strain_counter = 0
	empty_file_counter = 0

	summary_dict = defaultdict(dict)

	print('Parsing microbecensus output files...')

	for file in file_list:

		MG_ID = file.rstrip().split('.')[0]
		avg_genome_size = 0
		bases = 0
		genome_equivalents = 0

		file_counter += 1

		# parse file 
		if verbose:
			print(f'	- starting {file}')

		with open(file, 'r') as F:

			line_count = len(open(file).readlines())
			line_counter = 0

			for line in F:

				line_counter += 1

				if line_counter == 11:
					avg_genome_size=line.rstrip().split("\t")[1]
					# avg_genome_size=line.rstrip().split(":")[1]
					print(f'	   line {line_counter} : average genome size = {avg_genome_size}')

				elif line_counter == 12:
					bases=line.rstrip().split("\t")[1]
					# bases=line.rstrip().split(":")[1]
					print(f'	   line {line_counter} : bases = {bases}')

				elif line_counter == 13:
					genome_equivalents=line.rstrip().split("\t")[1]
					# genome_equivalents=line.rstrip().split(":")[1]
					print(f'	   line {line_counter} : genome equivalents = {genome_equivalents}')

				else:
					next


			summary_dict[f"{MG_ID}"] = {'MG_ID' : MG_ID,
										'avg_genome_size' : avg_genome_size,
										'bases' : bases,
										'genome_equivalents' : genome_equivalents}

	if verbose:
		print("")        
		print("summary dictionary:")
		print(summary_dict)

	print('')
	print(f' - parsed {file_counter} of {len(file_list)} microbe census output (*.census) files.')
	print(f' - wrote {len(summary_dict)} strain results to summary dictionary.')
	print('')

	return summary_dict


def dict_writer(summary_dict, output, verbose):
	
	print('-----------------------------------------------------------------------')
	print(f'Writing output with {len(summary_dict)} microbe census summaries...')

	# write output files 

	fieldnames = ('MG_ID', 'avg_genome_size', 'bases', 'genome_equivalents')

	outfile=f"{output}"
	with open(outfile, 'w', newline='') as f:
		writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
		writer.writeheader()
		for k in summary_dict:
			writer.writerow({field: summary_dict[k].get(field) or k for field in fieldnames})

	print(f' - output written to {output}')

def main():
	# configure argparse arguments & pass to functions.
	parser = argparse.ArgumentParser(
		description=__doc__,
		formatter_class = argparse.RawDescriptionHelpFormatter
		)
	parser.add_argument(
		'-i', '--input_directory',
		help = 'Directory containting *.census files.',
		metavar = '',
		type=str,
		required=True
		)	
	parser.add_argument(
		'-o', '--output',
		help = 'Please specify output file name (& optional prefix) with full path.',
		metavar = '',
		type=str,
		required=True
		)
	parser.add_argument(
		'-v', '--verbose',
		help = 'Toggle volume of printed output.',
		action='store_true'
		)
	args=vars(parser.parse_args())

	if args['verbose']:
		print('')
		print(f"Summary directory: {args['input_directory']}")
		print('')

	summary_dict = parse_census(args['input_directory'], args['verbose'])

	dict_writer(summary_dict, args['output'], args['verbose'])

	print('Looks like everything completed!')

if __name__ == "__main__":
	main()
