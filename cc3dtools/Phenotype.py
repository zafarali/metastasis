## JUNE 2015

from Genome import Genome, Mutation
from collections import Counter

class PhenotypeEvaluator( object ) :

	def __init__ ( self , phenotypes ) :
		"""
			PhenotypeEvaluator
				This class allows us to create a phenotype evaluator
				that accepts a dict of phenotypes on initialization.
				It can then evaluate individual Genome objects and
				return the phenotypes that it exhibits

			phenotypes / dictionary [mandatory]
				must contain entries akin to { 'phenotype-name': (r_low, r_high) , ... }
				where r_high > r_low, determine the region within which the loci 
				will be matched against to return 'phenotype-name'
		"""
		for k , v in phenotypes.items():
			assert type( k ) is str , 'phenotype keys must be strings'
			assert v[1] > v[0] , 'upper bound of ' + k + ' must be greater than the lower bound'
			assert type( v[1] ) is int and type( v[0] ) is int, ' (!) recent change means bounds need to be in ints now: https://github.com/zafarali/metastasis/issues/17'

		self.phenotypes = phenotypes
		

	def evaluate ( self , genome ) :
		"""
			evaluates phenotypes against a given genome
			@params:
				genome / Genome [mandator]
					genome against which we are going to check loci for phenotypic matches
			@return:
				list
					a list of phenotypes for which this genome has been matched
		"""

		assert isinstance( genome , Genome ), 'genome supplied must be of type cc3dtools.Genome!'
		loci = genome.get_mutated_loci()
		matched_phenotypes = []
		phenotypes = self.phenotypes.items()

		for locus in loci:
			for phenotype, region in phenotypes:
				# check if the locus is in the region
				# 'locus.locus' to get the float value of that mutation rather 
				# than an object!
				if locus.locus > region[0] and locus.locus < region[1]:
					matched_phenotypes.append( phenotype )
		return Counter( matched_phenotypes )


class Phenotype( object ):
	def __init__ ( self , phenotypes ):
		"""
			holds the phenotypes related to an individual

			phenotypes / dictionary [mandatory]
				must contain entries akin to { 'phenotype-name': (r_low, r_high) , ... }
				where r_high > r_low, determine the region within which the loci 
				will be matched against to return 'phenotype-name'
		"""
		self.counts = {}
		for k , v in phenotypes.items():
			assert type( k ) is str , 'phenotype keys must be strings'
			
			self.counts[ k ] = 0

			assert v[1] > v[0] , 'upper bound of ' + k + ' must be greater than the lower bound'

		self.phenotypes = phenotypes
		

	def evaluate ( self, mutation ) :
		"""
			evaluates the phenotype of a given Mutation object or
			mutation loci (int)

			mutation / Mutation or int [mandatory]
			mutation to evaluate
		"""
		if isinstance( mutation , Mutation ):
			mutation = mutation.to_int()

		assert type( mutation ) is int , 'mutation must work out to a int or have an internal representaiton of int'


		for phenotype, region in self.phenotypes.items():

			if mutation > region[0] and mutation < region[1]:
				# print phenotype
				self.counts[ phenotype ] += 1

	def get_counts ( self ):
		return self.counts

