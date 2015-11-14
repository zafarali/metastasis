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
		

	def replicate(self, name = ''):
		"""
			Replicates the chromosome and returns a new Chromosome object
			that is identical to the current one.
			@params:
				name: name of the new chromosome, defaults to ''
			@return:
				returns a new Chromosome() that is identical to this.
		"""

		replicated_chromosome = Chromosome ( mean_mutations = self.mean_mutations , chromosome_order = self.chromosome_order , name = name)

		replicated_chromosome.mutated_loci = replicated_chromosome.mutated_loci.union( self.mutated_loci )

		return replicated_chromosome

	def get_mutated_loci(self):
		"""
			Returns a set() of loci that have been mutated
		"""

		return self.mutated_loci

	def is_mutated(self, locus):
		"""
			Checks if a locus is mutated
			@params:
				locus: the position of the gene
		"""

		if isinstance( locus, Mutation):
			return locus in self.mutated_loci

		return Mutation(locus) in self.mutated_loci


class Genome(object):
	def __init__(self, mean_mutations=0, chromosome_order=15, ploidy=2, name='GenericGenome'):
		self.chromosomes = [ Chromosome(mean_mutations=mean_mutations, chromosome_order=chromosome_order,name=name+'_'+str(i)) for i in xrange(ploidy) ]
		pass
	def mutate(self):
		"""
			Mutates all chromosomes in this genome
		"""
		for chromosome in self.chromosomes:
			chromosome.mutate()

	def replicate(self):

		replicated_genome = Genome(mean_mutations=self.mean_mutations, chromosome_order=self.chromosome_order, ploidy=self.ploidy, name=self.name+'r' )

		replicated_chromosomes = []

		for chromosome in self.chromosomes:
			replicated_chromosomes.append( chromosome.replicate() )

		replicated_genome.chromosomes = replicated_chromosomes

		return replicated_genome


class Phenotype(object):
	def __init__(self, phenotype_template):
		pass

class pDivisionFunction(object):
	@staticmethod
	def constant(p):
		def const_fn(*args,**kwargs):
			return p
		return const_fn

class Cell(object):
	def __init__(self, mean_mutations=0, chromosome_order=15, ploidy=2, name='GenericCell', cell_type=1, phenotype_specs=False, p_division_function=pDivisionFunction.constant(0.5)):
		self.cell_type = 1
		self.genome = Genome(mean_mutations=mean_mutations, chromosome_order=chromosome_order, name=name+'_g', ploidy=ploidy)
		self.phenotype = Phenotype( {} if not phenotype_specs else phenotype_specs )
		self.p_division_function = p_division_function
		self.attributes = {}

	def p_division(self, **kwargs):
		"""
			Returns the probability of this cell dividing
		"""
		return self.p_division_function(self.attributes, **kwargs)
		# pass

class SelectionDistribution(object):
	@staticmethod
	def equal(cell_array):
		raw_dist = np.array([ cell.p_division() for cell in cell_array ])
		return raw_dist / float(np.sum(raw_dist))
	@staticmethod
	def cancer_only(cell_array):
		raw_dist = np.array( [ cell.p_division() if cell.is_cancer() else 0 for cell in cell_array ] )
		return raw_dist / float( np.sum(raw_dist) )

class Simulator(object):
	def __init__(self):
		pass

	def run(self, time_steps = 100):
		"""
			Runs the model for time_steps
		"""
		for i in xrange(time_steps):
			self.step()
		

	def step(self):
		print('Took Step')
		pass