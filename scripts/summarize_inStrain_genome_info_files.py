#!/usr/bin/env python

'''Parses inStrain genome coverage/breadth results to 
calculate summary statistics for each genome/MAG.

As input, this script requires:
	- path to directory of linked inStrain *genome_info.tsv files

The output produced is:
	- a tsv summary of coverage breadth and depth
	for the input *genome_info.tsv files

'''

import argparse, glob, os, csv
import pandas as pd
from statistics import mean, median
from collections import defaultdict

# define parse_genome_info function
def parse_genome_info(input_directory, verbose):
	# generate dictionary of genome coverage depth & breadth values

	os.chdir(input_directory)					# change to input directory
	file_list = glob.glob("*genome_info.tsv")			# generate list of *.census files

	file_counter = 0
	strain_counter = 0
	empty_file_counter = 0

	summary_dict = defaultdict(dict)

	print('Parsing inStrain genome_info.tsv files...')

	for file in file_list:

		ID = file.rstrip().split('-')[1]
		MG_visit = file.rstrip().split('-')[2]
		MG_ID = f"{ID}-{MG_visit}"
		breadth = 0
		depth = 0

		file_counter += 1

		# parse file 
		if verbose:
			print(f'	- starting {file}')

		with open(file, 'r') as F:
			next(F)

			strain_counter = 0

			for line in F:
				genome = line.rstrip().split('\t')[0]
				depth  = line.rstrip().split('\t')[1]
				breadth  = line.rstrip().split('\t')[2]

				######### will need revisiting
				######### because of variable field numbers output by inStrain

				# conANI_reference = line.rstrip().split('\t')[13]
				# popANI_reference = line.rstrip().split('\t')[14]

				# consensus_divergent_sites = line.rstrip().split('\t')[21]
				# population_divergent_sites = line.rstrip().split('\t')[22]
				# divergent_site_count = line.rstrip().split('\t')[29]


				if verbose:
					print(f'	   MG_ID {MG_ID} : genome = {genome} | breadth = {breadth}   depth = {depth}')


				summary_dict[f"{MG_ID}_{strain_counter}"] = {'MG_ID' : MG_ID,
											'genome' : genome,
											'depth' : depth,
											'breadth' : breadth}

				strain_counter += 1

				######### will need revisiting
				######### because of variable field numbers output by inStrain
			# summary_dict[f"{MG_ID}"] = {'MG_ID' : MG_ID,
			# 							'genome' : genome,
			# 							'depth' : depth,
			# 							'breadth' : breadth,
			# 							'conANI_reference' : conANI_reference,
			# 							'popANI_reference' : popANI_reference,
			# 							'consensus_divergent_sites' : consensus_divergent_sites,
			# 							'population_divergent_sites' : population_divergent_sites,
			# 							'divergent_site_count' : divergent_site_count}


	if verbose:
		print("")
		print("summary dictionary:")
		print(summary_dict)

	print('')


	print('')
	print(f' - parsed {file_counter} of {len(file_list)} inStrain output (*genome_info.tsv) files.')
	print(f' - wrote {len(summary_dict)} strain results to summary dictionary.')
	print('')

	return summary_dict


def dict_writer(summary_dict, output, verbose):
	
	print('-----------------------------------------------------------------------')
	print(f'Writing output with {len(summary_dict)} inStrain *genome_info.tsv summaries...')

	# write output files 

	fieldnames = ('MG_ID', 'genome', 'depth', 'breadth')


				######### will need revisiting
				######### because of variable field numbers output by inStrain
	# fieldnames = ('MG_ID', 'genome', 'depth', 'breadth',
	# 				'conANI_reference', 'popANI_reference',
	# 				'consensus_divergent_sites', 'population_divergent_sites',
	# 				'divergent_site_count')


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
		help = 'Directory containting *genome_info.tsv files.',
		metavar = '',
		type=str,
		required=True
		)	
	parser.add_argument(
		'-o', '--output',
		help = 'Please specify output file name with full path.',
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

	summary_dict = parse_genome_info(args['input_directory'], args['verbose'])

	dict_writer(summary_dict, args['output'], args['verbose'])

	print('Looks like everything completed!')

if __name__ == "__main__":
	main()
