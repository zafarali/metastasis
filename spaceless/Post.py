from collections import Counter
import numpy as np


def frequency_analyze( genomes, subsample=0 ):
	"""
		Returns the frequency distribution of mutations in the genomes
		@params:
			genomes: the genomes we wish to analyze
			subsample: the random subsample size of genomes we wish to pick
	"""
	counter = Counter()
	# create a random sample
	
	subsample = min( len( genomes ), subsample ) if subsample > 0 else len( genomes )

	idx = range( len( genomes ) )

	selected_idx = np.random.choice( idx, size=subsample, replace=False )

	

	for index in selected_idx:
		mutated_loci = genomes[ index ].get_mutated_loci()
		counter.update( mutated_loci )

	return Counter( [ v for _, v in counter.most_common() ] ), subsample


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
	return dict(Epi = Epi, S=S, SH=SH, D=D)

	




