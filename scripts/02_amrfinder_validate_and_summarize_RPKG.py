#!/usr/bin/env python

'''Validate and summarize RPKG of AMRFinder-detected genes 
for plots & analysis.

This script takes the following inputs:
- directory containing filtered AMRFinder tsv files
	(output from step 00)
- directory containing ${uniqueID}_gene_info.tsv files
	(output from inStrain https://instrain.readthedocs.io/)

With intermediate validation steps (option -V):

- all genes input in AMRFinder tables are tested against all genes 
in the *gene_info tsv files. If there are any genes that are 
not in the submitted gene_info tsv files, these are optionally 
output as: genes_to_validate.tsv

- all duplicated gene RPKG values are summed by sample.  
Input contigs/scaffolds hosting the detected genes are listed with
the summed (deduplicated) RPKG values and gene sequence name and
output as: deduplicated_RPKG.tsv

Genes that have RPKG values are returned with following output:
- specified output file & path containting a single tsv file with 
length / genome equivalent normalized AMR gene abundance (RPKG)

'''

import argparse, csv, glob, os
from collections import defaultdict
import pandas as pd
import glob, os


def parse_amrfinder_tsvs(input_directory, verbose, short_names, gene_metadata):
	# generate dictionary of genes & AMR sequence names

	# initialize dictionary & list objects
	os.chdir(input_directory)			# change to input directory
	file_list = glob.glob("*.tsv")		# generate list of *.tsv files
	AMR_dict = defaultdict(defaultdict(list).copy)		# initialize AMR_dict
	AMR_metadata_dict = defaultdict(defaultdict(list).copy)		# initialize AMR_metadata_dict
	MG_IDs = []

	print('Parsing AMRFinder tables...')


	for file in file_list:

		# mg_id = file.rstrip().split('.')[0]
		mg_id = file.rstrip().split('_')[0]
		if mg_id not in MG_IDs:
			MG_IDs.append(mg_id)

		with open(file, 'r') as F:
			for line in F:
				X = line.rstrip().split('\t')
				gene = X[0]				# gene ID / scaffold

				# pull long vs short gene names
				if short_names:
					sequence_name = X[1] 	# sequence name

				else:
					sequence_name = X[2] 	# sequence name

			# parse AMRFinder gene annotation metadata
				short_name = X[1]
				long_name = X[2]
				scope = X[3]
				element_type = X[4]
				element_subtype = X[5]
				gene_class = X[6]
				gene_subclass = X[7]
				method_type = X[8]

				if gene not in AMR_dict.keys():
					AMR_dict[gene]=sequence_name

				if gene_metadata:
					if long_name not in AMR_metadata_dict.keys():
						AMR_metadata_dict[long_name]={'short_name' : short_name,
														'long_name' : long_name,
														'scope' : scope,
														'element_type' : element_type,
														'element_subtype' : element_subtype,
														'gene_class' : gene_class,
														'gene_subclass' : gene_subclass,
														'method_type' : method_type}


	return AMR_dict, MG_IDs, AMR_metadata_dict


def parse_genome_equivalents(census_summary, verbose, AMR_dict):
	# generate dictionary of metagenome genome equivalents

	# initialize dictionary & list objects
	ge_dict = defaultdict(defaultdict(list).copy)		# initialize ge_dict

	print('Parsing MicrobeCensus summary table...')


	with open(census_summary, 'r') as F:
		for line in F:
			X = line.rstrip().split('\t')
			MG_ID = X[0]				# MG_ID
			# avg_genome_size = X[1]	# metagenome average genome size
			# bases = X[2]				# bases per metagenome
			genome_equivalents = X[3]	# genome equivalents per metagenome

			ge_dict[MG_ID]=genome_equivalents

	return ge_dict


def parse_coverage_tsvs(input_directory, verbose, AMR_dict):
	# generate dictionary of prodigal gene coverage & length

	# initialize dictionary & list objects
	os.chdir(input_directory)					# change to input directory
	file_list = glob.glob("*_gene_info.tsv")	# generate *.tsv list
	RPKG_dict = defaultdict(defaultdict(list).copy)		# initialize RPKG_dict
	length_dict = defaultdict(defaultdict(list).copy)		# initialize length_dict

	print('Parsing inStrain coverage tables...')


	for file in file_list:

		with open(file, 'r') as F:
			next(F)

			for line in F:
				X = line.rstrip().split('\t')

				# try to skip lines with empty length values
				if X[2] != '':

					# print(f"X[2]: {X[2]}")

					gene = X[1]			# prodigal gene name / scaffold

					coverage = float(X[3]) 	# coverage depth value
					# nucl_diversity = X[6]	# nucleic acid diversity
					# SNV_count = X[13]	# SNV counts for gene

					length = float(X[2])	# gene length

					if gene in AMR_dict.keys():
						RPKG_dict[gene]=coverage
						length_dict[gene]=length

					# print(f"gene: {gene}	length: {length}")

	return RPKG_dict, length_dict


def validate(RPKG_dict, AMR_dict, AMR_input_directory, verbose):
	# test differences in genes with RPKG values vs those without

	# initialize dictionary & list objects
	validate_dict = defaultdict(defaultdict(list).copy)		# initialize validate_dict
	validate_detail_dict = defaultdict(list)		# initialize validate_dict
	total_genes = 0
	valid_gene = 0
	method_hmm = 0

	mg_to_validate = []
	genes_to_validate = []
	file_list = []
	RPKG_genes = RPKG_dict.keys()

	print('Validating results...')


	for gene, sequence_name in AMR_dict.items():
		total_genes += 1

		if gene not in RPKG_dict.keys():
			# print(f"	  gene to validate: {gene}")
			# print(f"	  sequence_name to validate: {sequence_name}")
			validate_dict[gene]=sequence_name
		else:
			valid_gene += 1

	print('')
	print(f"   Number of genes checked: {total_genes}")
	print(f"   Number of valid genes: {valid_gene} ({round((valid_gene / total_genes * 100), 2)} %)")
	print(f"   Number of genes needing validation: {total_genes - valid_gene} ({round(((total_genes-valid_gene) / total_genes * 100), 2)} %)")
	print('')

	for gene in validate_dict.keys():
		mg_gene_id = gene.rstrip().split('_')[0]
		if mg_gene_id not in mg_to_validate:
			mg_to_validate.append(mg_gene_id)

	if verbose:
		print('')
		print('   Pulling list of metagenomes to validate...')
#		print(mg_to_validate)
		print('')

	os.chdir(AMR_input_directory)		# change to AMRFinder tsv input directory
	
	for mg in mg_to_validate:
		if glob.glob(f"{mg}*.tsv") not in file_list:
			file_list.append(glob.glob(f"{mg}*.tsv"))

	if verbose:
		print('   AMRFinder files to validate:')
		print(file_list)
		print('')

	for mg_file in file_list:
		with open(f"{AMR_input_directory}/{mg_file[0]}", 'r') as file:
			for line in file:
				X = line.rstrip().split('\t')
				protein_id = X[0]		# protein id (with contig)
				gene_symbol = X[1]		# gene symbol
				sequence_name = X[2] 	# sequence name
				scope = X[3]			# core / plus
				element_type = X[4]		# AMR, STRESS, etc.
				element_subtype = X[5]	# AMR, Metal, Acid, etc.
				amr_class = X[6]		# e.g. glycopeptide, aminoglycoside, etc.
				amr_subclass = X[7]		# e.g. vancomycin, streptomycin, etc.
				method_type = X[8]			# e.g. PARTIALP, HMM, EXACTP, etc.
				target_l = X[9]			# target length
				ref_l = X[10]			# ref length
				cov_pct = X[11]			# coverage percent [breadth]
				pct_ident = X[12]		# percent identity to reference
				align_l = X[13]			# alignment length
				accession = X[14]		# NCBI accession for closest sequence
				ref_name = X[15]		# name of closest sequence
				HMM_id = X[16]			# id of closest HMM
				HMM_desc = X[17]		# description of closest HMM


				# if protein_id in validate_dict.keys():
				# 	validate_detail_dict[protein_id]=line

				if protein_id in validate_dict.keys():
					validate_detail_dict[protein_id]=[protein_id.rstrip().split('_')[0], 
					gene_symbol, sequence_name]

	return validate_dict, validate_detail_dict


def generate_RPKG_matrix(AMR_dict, RPKG_dict, length_dict, ge_dict, MG_IDs, verbose, raw):
	print('Generating RPKG matrix...')

	# initialize dict, lists
	# AMR_dict - 	key: gene/scaffold name 	value: sequence name (unique)
	# RPKG_dict - 	key: gene/scaffold name  	value: RPKG value
	# length_dict -	key: gene/scaffold name 	value: gene length
	# ge_dict -		key: MG_ID 					value: genome equivalents

	AMR_scaffold_RPKG = defaultdict(list)
	scaf_mgid_name_rpkg = defaultdict(int)
	scaf_mgid_name_raw = defaultdict(int)
	matrix = defaultdict(list)
	raw_matrix = defaultdict(list)
	mg_genes = defaultdict(list)
	raw_mg_genes = defaultdict(list)
	dedup_dict = defaultdict(list)
	raw_dedup_dict = defaultdict(list)

	unique_gene_list = []

	# generate list of unique AMR genes
	for scaffold, gene_name in AMR_dict.items():
		if gene_name not in unique_gene_list:
			unique_gene_list.append(gene_name)
		# mg_id=scaffold.rstrip().split('_')[0]
		# print(f' gene_name: {gene_name},   MG_ID: {mg_id}')

	# output counts if verbose
	if verbose:
		print(f"{len(MG_IDs)} MG IDs available | {len(unique_gene_list)} unique genes detected")

	# merge AMR_dict & RPKG_dictionaries as scaffold : gene_name, RPKG
	for d in (AMR_dict, RPKG_dict, length_dict):
		for key, value in d.items():
			AMR_scaffold_RPKG[key].append(value)

	# split into structured dictionary - scaffold : gene name, RPKG,  MG ID
	for scaffold, name_RPKG in AMR_scaffold_RPKG.items():
		mg_id = scaffold.rstrip().split('_')[0]
		g_eqs = float(ge_dict[mg_id])

		# print(f'genome_equivalents: {g_eqs}')

		if len(name_RPKG) > 1:

			# print(f'name_RPKG[1]: {name_RPKG[1]}')
			# print(f'name_RPKG[2]: {name_RPKG[2]}')
			g_length = round(name_RPKG[2])

			scaf_mgid_name_rpkg[scaffold] = { 'name' : name_RPKG[0],
											  'RPKG' : float(name_RPKG[1])/g_length/g_eqs,
											  # RPKG calculated as mapped reads / gene length / genome equiv.
											  # see  https://github.com/snayfach/MicrobeCensus

											  'mg_id' : mg_id
											}

			if raw:
				scaf_mgid_name_raw[scaffold] = { 'name' : name_RPKG[0],
												  'raw_coverage' : name_RPKG[1],
												  'mg_id' : mg_id
												}

			print(f"name_RPKG > 1 | name: {name_RPKG[0]}  RPKG: {float(name_RPKG[1])/g_length/g_eqs}")
			print(f"	- coverage: {name_RPKG[1]}	gene length: {g_length}	genome eqs: {g_eqs}")


		else:
			scaf_mgid_name_rpkg[scaffold] = { 'name' : name_RPKG[0],
											  'RPKG' : 0,
											  'mg_id' : scaffold.rstrip().split('_')[0]
											}
			# print(f"name_RPKG <= 1 | name: {name_RPKG[0]}  RPKG: {scaf_mgid_name_rpkg[scaffold]['RPKG']}")


			if raw:
				scaf_mgid_name_raw[scaffold] = { 'name' : name_RPKG[0],
												  'raw_coverage' : 0,
												  'mg_id' : mg_id
												}


	# tally & deduplicate AMR gene RPKG values by sample
	for MG in MG_IDs:
		dedup_list = defaultdict(list)
		raw_dedup_list = defaultdict(list)
		dedup_count = 0

		if verbose:
			print('')
			print(f"Evaluating gene duplicates in {MG}...")

		for scaffold, values in scaf_mgid_name_rpkg.items():
			# if verbose:
			# 	print(f"scaffold: {scaffold}  mg_id: {values['mg_id']} RPKG: {values['RPKG']}")

			if MG == values['mg_id'] and values['name'] not in dedup_list.keys():
				dedup_list[values['name']]={'mg_id' : values['mg_id'],
											'RPKG': values['RPKG'], 
											'count' : 1, 
											'scaffolds' : [scaffold],
											'gene_name' : values['name']}
				# print(dedup_list[values['name']])

			elif MG == values['mg_id'] and values['name'] in dedup_list.keys():
				dedup_list[values['name']]['RPKG']  = float(values['RPKG']) + float(dedup_list[values['name']]['RPKG'])
				dedup_list[values['name']]['count'] = int(dedup_list[values['name']]['count']) +1
				dedup_list[values['name']]['scaffolds'].append(scaffold)
				dedup_count += 1

			if raw:
				if MG == values['mg_id'] and values['name'] not in raw_dedup_list.keys():
					# print(f" raw_coverage: {scaf_mgid_name_raw[scaffold]['raw_coverage']}")
					raw_dedup_list[values['name']]={'mg_id' : values['mg_id'],
												'raw_coverage': scaf_mgid_name_raw[scaffold]['raw_coverage'],
												'count' : 1, 
												'scaffolds' : [scaffold],
												'gene_name' : values['name']}

				elif MG == values['mg_id'] and values['name'] in raw_dedup_list.keys():
					raw_dedup_list[values['name']]['raw_coverage']  = scaf_mgid_name_raw[scaffold]['raw_coverage'] + raw_dedup_list[values['name']]['raw_coverage']
					raw_dedup_list[values['name']]['count'] = int(dedup_list[values['name']]['count']) +1
					raw_dedup_list[values['name']]['scaffolds'].append(scaffold)
					dedup_count += 1


		if verbose:
			#print(dedup_list)
			print(f"Number of deduplicated genes: {dedup_count}")	
			for gene_name, value in dedup_list.items():
				if value['count'] > 1:
					print('')
					print(f"Duplicate gene is: '{gene_name}' count: {dedup_list[gene_name]['count']}")
					print(f"	 tallied RPKG: {dedup_list[gene_name]['RPKG']}")
					print(f"	 tallied raw coverage: {raw_dedup_list[gene_name]['raw_coverage']}")
					for scaffold, o_values in scaf_mgid_name_rpkg.items():
						if o_values['name'] == gene_name and o_values['mg_id'] == MG:
							print(F"		- original name '{o_values['name']}' | original RPKG {o_values['RPKG']}")

		# add to dedup_dict to store all values across samples
		for gene_name, value in dedup_list.items():
			mg_genes[MG].append(gene_name)
			dedup_dict[f"{MG}_{gene_name}"]=[value['mg_id'], value['RPKG'], 
			value['count'], value['scaffolds'], value['gene_name']]

		if raw:
			for gene_name, value in raw_dedup_list.items():
				raw_mg_genes[MG].append(gene_name)
				raw_dedup_dict[f"{MG}_{gene_name}"]=[value['mg_id'], value['raw_coverage'], 
				value['count'], value['scaffolds'], value['gene_name']]


		#add MG gene RPKG values to matrix
		print('')
		print(f"Updating matrix with {MG} values...")
		for gene in unique_gene_list:
			if gene in mg_genes[MG]:
				matrix[MG].append(dedup_list[gene]['RPKG'])
				# if verbose:
				# 	print(f"	 {gene} (RPKG {dedup_list[gene]['RPKG']}) added to matrix")
			else:
				matrix[MG].append(0)
				# if verbose:
				# 	print(f"	 {gene} (RPKG 0, ND) added to matrix")		

		if raw:
			print(f"Updating matrix with {MG} raw coverage values...")
			for gene in unique_gene_list:
				if gene in raw_mg_genes[MG]:
					raw_matrix[MG].append(raw_dedup_list[gene]['raw_coverage'])
					# if verbose:
					# 	print(f"	 {gene} (RPKG {dedup_list[gene]['RPKG']}) added to matrix")
				else:
					raw_matrix[MG].append(0)
					# if verbose:
					# 	print(f"	 {gene} (RPKG 0, ND) added to matrix")		


	# if verbose:
	#	 print('')
	#	 print("MG gene list:")
	#	 print(f"   length: {len(mg_genes)}")
	#	 print(mg_genes)
	#	 print('')
	#	 print("dedup_dict:")
	#	 print(f"   length: {len(dedup_dict)}")
	#	 print(dedup_dict)
	#	 print('')
	#	 print("matrix:")
	#	 print(f"   length: {len(matrix)}")
	#	 print(matrix)	


	print('')
	print('===========================================')
	print('Converting completed matrix to dataframe...')

	RPKG_matrix = pd.DataFrame(matrix, index = unique_gene_list)
	RPKG_matrix.sort_index(inplace=True)
	RPKG_matrix.sort_index(axis=1, inplace=True)	

	if raw:
		print('Converting completed raw coverage matrix to dataframe...')

		R_matrix = pd.DataFrame(raw_matrix, index = unique_gene_list)
		R_matrix.sort_index(inplace=True)
		R_matrix.sort_index(axis=1, inplace=True)	


	print(f"{len(MG_IDs)} MGs evaluated with {len(unique_gene_list)} unique genes.")
	print('===========================================')
	
	return RPKG_matrix, unique_gene_list, dedup_dict, R_matrix, raw_dedup_dict

def main():
	# configure argparse arguments & pass to functions.
	parser = argparse.ArgumentParser(
		description=__doc__,
		formatter_class = argparse.RawDescriptionHelpFormatter
		)
	parser.add_argument(
		'-a', '--amrfinder_tsv_path',
		help = 'Please specify directory path containing filtered AMRFinder tsv files.',
		#metavar = '',
		type=str,
		required=True
		)
	parser.add_argument(
		'-i', '--inStrain_tsv_path',
		help = 'Please specify directory path containing inStrain gene_info tsv path.',
		#metavar = '',
		type=str,
		required=True
		)
	parser.add_argument(
		'-m', '--microbe_census_summary',
		help = 'Please specify microbe census summary file with full path.',
		#metavar = '',
		type=str,
		required=True
		)
	parser.add_argument(
		'-o', '--output',
		help = 'Please specify output file path (& optional prefix).',
		#metavar = '',
		type=str,
		required=True
		)
	parser.add_argument(
		'-s', '--short_names',
		help = 'Option to summarize AMR genes by short names.',
		action='store_true'
		),
	parser.add_argument(
		'-r', '--raw_coverage',
		help = 'Option to write raw (not-normalized) mean gene coverage tsv matrix.',
		action='store_true'
		),
	parser.add_argument(
		'-V', '--validate',
		help = 'Write genes_to_validate.tsv and deduplicated.tsv.',
		action='store_true'
		)
	parser.add_argument(
		'-g', '--gene_metadata',
		help = 'Write gene_metadata.tsv.',
		action='store_true'
		)	
	parser.add_argument(
		'-v', '--verbose',
		help = 'Toggle volume of printed output.',
		action='store_true'
		)
	args=vars(parser.parse_args())

	if args['verbose']:
		print('')
		print('AMRFinder tsv directory:', args['amrfinder_tsv_path'])
		print('inStrain tsv directory:', args['inStrain_tsv_path'])
		print('MicrobeCensus summary:', args['microbe_census_summary'])
		print('')

	AMR_dict, MG_IDs, AMR_metadata_dict= parse_amrfinder_tsvs(args['amrfinder_tsv_path'], args['verbose'], args['short_names'], args['gene_metadata'])
	ge_dict= parse_genome_equivalents(args['microbe_census_summary'], args['verbose'], AMR_dict)
	RPKG_dict, length_dict= parse_coverage_tsvs(args['inStrain_tsv_path'], args['verbose'], AMR_dict)
	validate_dict, validate_detail_dict= validate(RPKG_dict, AMR_dict, args['amrfinder_tsv_path'], args['verbose'])	
	RPKG_matrix, unique_gene_list, dedup_dict, raw_matrix, raw_dedup_dict= generate_RPKG_matrix(AMR_dict, RPKG_dict, length_dict, ge_dict, MG_IDs, args['verbose'], args['raw_coverage'])

# option to write validation tsv files
	if args['validate']:
		print('')
		print(f"Writing validation files...")
		print(f"Validation output path: {args['output']}")
		def list_to_tsv(path, out_file_name, dictionary,
			col_header):
			outfile=f"{path}/{out_file_name}"
			df = pd.DataFrame.from_dict(dictionary,
				orient="index", columns=col_header)
			df.to_csv(outfile,  sep='\t')
			print(f"   ...{out_file_name} complete!")

		print('')
		print(f"   - validate_detail_dict...")

		col_header=['mg_id', 
					'gene_symbol', 
					'sequence_name']
		list_to_tsv(args['output'], 'genes_to_validate.tsv',
			validate_detail_dict, col_header)

		print('')
		print(f"   - dedup_dict...")

		col_header=['mg_id', 
			'RPKG', 'count', 'scaffolds', 'gene_name']
		list_to_tsv(args['output'], 'deduplicated_RPKG.tsv',
			dedup_dict, col_header)


		# validate_detail_dict.to_csv(f"{args['output']}/genes_to_validate.tsv", sep='\t')
		# dedup_dict.to_csv(f"{args['output']}/deduplicated.tsv", sep='\t')		

# option to write gene metadata tsv file
	if args['gene_metadata']:
		print('')
		print(f"Writing gene metadata files...")
		print(f"Gene metadata output path: {args['output']}")

		fieldnames = ('short_name', 'long_name', 'scope', 'element_type', 
						'element_subtype', 'gene_class', 'gene_subclass',
						'method_type')

		outfile=f"{args['output']}/gene_metadata.tsv"
		with open(outfile, 'w', newline='') as f:
			writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
			writer.writeheader()
			for k in AMR_metadata_dict:
				writer.writerow({field: AMR_metadata_dict[k].get(field) or k for field in fieldnames})

		if args['verbose']:
			print('')
			print(f"Gene metadata parsed for {len(AMR_metadata_dict)} genes")
			print(f"Gene metadata written to {outfile}")


# option to write raw mean coverage tsv file
	if args['raw_coverage']:
		print('')
		print(f"Writing raw coverage tsv files...")
		raw_matrix.to_csv(f"{args['output']}/raw_coverage_matrix.tsv", sep='\t')


# write output tsv file
	RPKG_matrix.to_csv(f"{args['output']}/RPKG_matrix.tsv", sep='\t')

	print('===========================================')
	print('')
	print('Looks like everything completed!')
	print(f"Output files written to: {args['output']}")

	if args['verbose']:
		#print('')
		#print('AMR dictionary:')
		#print(AMR_dict)
		# print('')
		# print('RPKG dictionary:')
		# print(RPKG_dict)
		# print('')
		# print('Unique gene list:')
		# print(unique_gene_list)
		# print(len(unique_gene_list), 'unique genes detected.')
		print('')
		print(len(AMR_dict),'Genes with AMRFinder hits were identified.')
		print(len(RPKG_dict),'Genes with inStrain coverage values were identified.')

	#print('Complete. Output written to:')
	#print(args['output'])

if __name__ == "__main__":
	main()
