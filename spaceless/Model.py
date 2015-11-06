import random
import numpy as np

class Mutation(object):
	def __init__(self, locus):
		self.locus = int( locus )

	def __repr__(self):
		"""
			mutation objects are represented as #LOCUS
		"""
		return '#'+str( self.locus )

	def __str__(self):
		"""
			when converted to a string, mutation objects are just represented by their loci
		"""
		return str( self.locus )

	def to_int(self):
		return self.locus

	def __cmp__(self, other):
		# print self.locus, other.locus
		if self.locus > other.locus:
			# print 'big'
			return 1
		elif self.locus < other.locus:
			# print 'small'
			return -1
		elif  self.locus == other.locus:
			# print 'equal'
			return 0

	def __hash__(self):
		return hash(self.locus)

class Chromosome(object):
	def __init__(self, mean_mutations=0, chromosome_order=15, name='GenericChromosome'):
		"""
			mutation_rate / int / 0
				the mean number of new mutations that will be produced
				every mutate() event.
			chromosome_order / int / 15
				the number of loci in the genome ( 10^chromosome_order )
			name / str / ''
				the name of the chromosome.
		"""

		self.mutated_loci = set()
		self.mean_mutations = mean_mutations
		self.chromosome_order = chromosome_order
		self.name = name

	def mutate(self):
		"""
			Introduces mutations into the chromosome
			@return:
				set() of loci which were mutated in this division
		"""
		number_of_new_mutations = np.random.poisson( self.mean_mutations )
		number_of_loci = 10**self.chromosome_order

		new_mutated_loci = set( map( Mutation , np.random.randint( number_of_loci , size = number_of_new_mutations ) ) )

		self.mutated_loci = self.mutated_loci.union( new_mutated_loci )

		return new_mutated_loci
		

	def replicate(self, name):
		"""
			Replicates the chromosome and returns a new Chromosome object
			that is identical to the current one.
		"""
		pass


class Genome(object):
	def __init__(self):
		pass

class Phenotype(object):
	def __init__(self):
		pass

class Cell(object):
	def __init__(self):
		pass

class Simulator(object):
	def __init__(self):

		# create the first cell
		self.cell_ids = [ Cell() ]
		self.genomes = { 
			0: Genome() 
		}
		self.phenotypes = { 
			0: Phenotype() 
		}

	def run(self, time_steps = 100):
		"""
			Runs the model for time_steps
		"""
		for i in xrange(time_steps):
			self.step()
		

	def step(self):
		print('Took Step')
		pass