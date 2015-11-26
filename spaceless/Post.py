from collections import Counter
import numpy as np
import random

def frequency_analyze( genomes, subsample = 0, threshold=0, cell_types=[] ):
	"""
		Returns the frequency distribution of mutations in the genomes
		@params:
			genomes: the genomes we wish to analyze
			subsample: the random subsample size of genomes we wish to pick
	"""
	counter = Counter()
	# create a random sample
	
	subsample = min( len( genomes ), subsample ) if subsample > 0 else len( genomes )


	iteration = 0
	while True:
		idx = range( len( genomes ) )
		selected_idx = np.random.choice( idx, size=subsample, replace=False )
		# assert threshold !=0 and cell_types

		if threshold == 0: break
		if iteration > 10: 
			print('Thresholding could not be met')
			break

		cell_types = np.array(cell_types)
		selected_cell_types = cell_types[selected_idx]
		cell_counts = Counter( selected_cell_types.tolist() )

		normal = cell_counts[1]
		cancer = cell_counts[2]

		proportion_cancer = float(cancer) / float(normal + cancer)
		iteration += 1
		if proportion_cancer > threshold: break


	for index in selected_idx:
		mutated_loci = genomes[ index ].get_mutated_loci()
		counter.update( mutated_loci )

	return Counter( [ v for _, v in counter.most_common() ] ), subsample

# def types_

def proportion_of_pairwise_differences ( allele_frequencies, number_of_genomes ):
	"""
		Calculates E(pi) = the expected proportion of pairwise
						differences 
		@params:
			allele_frequencies:
				the output from PostProcess.sample_analyze(sample)[:][1]
				(i.e not the distances, loop over each sample analyzed.)
		@return:
			E_of_pi : calculated according to the equation:
						## E(pi) = sum_{l=1}^M 2 * X_l ( 1 - X_l)
						## l = locus , M = number of loci
						## X_l = pairwise proportional difference at site l

	"""



	if len( allele_frequencies ):
		# 
		number_of_mutations, frequency_of_mutations = zip( *allele_frequencies.items() )
		# print 'frequency_of_mutations=',frequency_of_mutations
		number_of_mutations = np.array(number_of_mutations)  

		x = number_of_mutations / (2 * float(number_of_genomes) ) # proportion of cells in cluster mutated.
		# multiply by 2 to account for diploidy

		y = ( 2 * float(number_of_genomes) - number_of_mutations ) / ( 2 * float(number_of_genomes) - 1 )

		## E(pi) = sum_{l=1}^M 2 * X_l ( 1 - X_l)
		## l = locus , M = number of loci
		## X_l = pairwise proportional difference at site l
		individual_pis = 2 * x * y * np.array(frequency_of_mutations)
		# print 'individual pis=', individual_pis
		E_of_pi = np.sum( individual_pis ) 

		return E_of_pi 

	return 0

def number_of_segregating_sites ( allele_frequencies, *args, **kwargs ):
	"""
		Calculates S = number of segregating sites
		@params:
			allele_frequencies:
				the output from PostProcess.sample_analyze(sample)[:][1]
				(i.e not the distances, loop over each sample analyzed.)
			normalized:
				returns S / H_(n-1)
		@return:
			S : number of segregating sites

	"""

	if len( allele_frequencies ):
		# 
		number_of_mutations, frequency_of_mutations = zip( *allele_frequencies.items() )
		# print 'frequency_of_mutations=',frequency_of_mutations
		# print 'number_of_mutations=',number_of_mutations
		# print '',f
		number_of_mutations = np.array(number_of_mutations)  

		# multiply by 2 to account for diploidy

		S = np.sum( np.ones( len(number_of_mutations) - 1 ) * np.array(frequency_of_mutations[:-1]) ) 

		if kwargs.get('normalized', False):
			S = S / H(args[0]-1)
		return S

	return 0
	

def H(n):
	"""Returns an approximate value of n-th harmonic number.

	   http://en.wikipedia.org/wiki/Harmonic_number
	"""
	# Euler-Mascheroni constant
	gamma = 0.57721566490153286060651209008240243104215933593992
	return gamma + np.log(n) + 0.5/n - 1./(12*n**2) + 1./(120*n**4)

def get_stats(allele_frequencies, number_of_genomes):
	"""
		calculates statistics using allele_frequency distribution
		@returns:
			tuple containing:
				Epi
				S
				SH
				D
	"""
	S = number_of_segregating_sites(allele_frequencies, number_of_genomes, normalized=False)
	SH = S / H(number_of_genomes-1)
	Epi = proportion_of_pairwise_differences(allele_frequencies, number_of_genomes)
	D = Epi-SH
	return (S, SH, Epi, D)

	

def split_genomes(genomes, N, t=0):
	"""
		Generates a subsample of the genomes to 
		@params:
			genomes: genomes to shufffle and pick from:
				format: {1: [genome1, genome2, ...], 2: [genome3, genome4, ...] }

			N: the total number of genomes to return
			t: the thresholding level (% of cancer genomes in the sample)
	"""
	available_sample_size = len(genomes['cancer'])+len(genomes['normal'])

	if N > available_sample_size:
		return random.sample(genomes['cancer'] + genomes['normal'], available_sample_size), available_sample_size

	# number of cancer cells wanted in the final sample
	num_cancer_cells = int(N*t)


	# pick the final subsample size based on the minimum
	subsample_cancer = min( len(genomes['cancer']), num_cancer_cells )
	cancer_genomes = random.sample(genomes['cancer'],  subsample_cancer)

	# the number of normal cells need to make up the remaining cells
	num_normal_cells = N-subsample_cancer

	# pick the final subsample size based on the minimum.
	subsample_normal = min(len(genomes['normal']), num_normal_cells)
	normal_genomes = random.sample(genomes['normal'], subsample_normal)

	return normal_genomes + cancer_genomes, len(normal_genomes + cancer_genomes)


