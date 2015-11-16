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

	def mutate(self, number_of_new_mutations = False):
		"""
			Introduces mutations into the chromosome
			@params:
				number_of_new_mutations
					the number of mutations to create (by default it is poisson)
			@return:
				set() of loci which were mutated in this division
		"""
		if not number_of_new_mutations:
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
		self.mean_mutations = mean_mutations
		self.chromosome_order = chromosome_order
		self.ploidy = ploidy
		self.name = name
		pass
	def mutate(self, number_of_new_mutations=False):
		"""
			Mutates all chromosomes in this genome
		"""
		new_mutated_loci = set()
		for chromosome in self.chromosomes:
			new_mutated_loci  = new_mutated_loci.union( chromosome.mutate(number_of_new_mutations=number_of_new_mutations) )

		return new_mutated_loci

	def replicate(self, name=False):
		"""
			Returns a replicate of this genome.
		"""
		if not name:
			name = self.name
		replicated_genome = Genome(mean_mutations=self.mean_mutations, chromosome_order=self.chromosome_order, ploidy=self.ploidy, name=name )

		replicated_chromosomes = []

		for chromosome in self.chromosomes:
			replicated_chromosomes.append( chromosome.replicate() )

		replicated_genome.chromosomes = replicated_chromosomes

		return replicated_genome
	def get_mutated_loci(self):
		all_loci = set()
		for chromosome in self.chromosomes:
			all_loci = all_loci.union( chromosome.get_mutated_loci() )

		return all_loci

class Phenotype(object):
	def __init__(self, phenotype_template):
		"""
			holds the phenotypes related to an individual

			phenotype_template / dictionary
				must contain entries akin to { 'phenotype-name': (r_low, r_high) , ... }
				where r_high > r_low, determine the region within which the loci 
				will be matched against to return 'phenotype-name'
		"""
		self.counts = {}

		for k, v in phenotype_template.items():
			assert type(k) is str, 'phenotype keys  must be strings'

			self.counts[k] = 0

			assert v[1] > v[0], 'upper bound of ' + k + ' must be greater than lower bound '
		
		self.phenotypes = phenotype_template

	def _evaluate( self, mutation ):
		if isinstance( mutation , Mutation ):
			mutation = mutation.to_int()

		assert type( mutation ) is int , 'mutation must work out to a int or have an internal representaiton of int'


		for phenotype, region in self.phenotypes.items():

			if mutation > region[0] and mutation < region[1]:
				# print phenotype
				self.counts[ phenotype ] += 1

	def evaluate( self, mutations ):
		"""
			Evaluates a list of mutations to increase their counts in the phenotype
		"""


		if type(mutations) in [ list, set, tuple ]:
			for mutation in mutations:
				self._evaluate( mutation )
		else:
			self._evaluate(mutations)

	def get_counts ( self ):
		return self.counts

	def replicate ( self ):
		"""
			returns a copy of the current phenotype counter
		"""
		to_be_returned = Phenotype ( self.phenotypes )
		to_be_returned.counts = dict( self.get_counts().items() )

		return to_be_returned

class pDivisionFunction(object):
	"""
		Contains division functions
	"""
	@staticmethod
	def constant(p):
		def const_fn(*args,**kwargs):
			return p
		return const_fn

	@staticmethod
	def sigmoid(a=0.01,c=0):
		"""
			probability of division is proportional 
			to the number of mutations
		"""
		def sig_fn(*args, **kwargs):
			# let x = number of advantageous mutations 
			x = float(args[1].get_counts()['advantageous'])
			return 1. / ( 1 + np.exp( -( a*x - c ) ) )

		return sig_fn


# y_intercept = np.log( 1 / float( GLOBAL['_dV'] ) -1 )

# m = phenotypes[cell.id].get_counts()['advantageous']
# denominator = 1 + np.exp( -( 0.01*m - y_intercept ) )
# cell.targetVolume = min ( 1 / denominator + cell.targetVolume , GLOBAL['maxTargetVolume'] )



class Cell(object):
	def __init__(self, mean_mutations=0, chromosome_order=15, ploidy=2, name='GenericCell', cell_type=1, phenotypes=False, p_division_function=pDivisionFunction.constant(0.5), parent_name=None):
		self.cell_type = cell_type
		self.genome = Genome(mean_mutations=mean_mutations, chromosome_order=chromosome_order, name=name+'_g', ploidy=ploidy)
		self.phenotype = Phenotype( {} if not phenotypes else phenotypes )
		self.p_division_function = p_division_function
		self.attributes = {}
		self.name = name
		self.parent_name = parent_name
		self.generation = 0

	def p_division(self, **kwargs):
		"""
			Returns the probability of this cell dividing
		"""
		return self.p_division_function(self.attributes, self.phenotype, **kwargs)
		# pass

	def mitosis(self, name=False):
		"""
			Conducts a mitosis event and returns a new child cell
		"""
		if not name:
			name = self.name+'child'

		# create new cell
		new_cell = Cell(name=name, cell_type =self.cell_type, p_division_function=self.p_division_function)

		# mutate and replicate the genome
		new_mutations = self.genome.mutate()

		self.phenotype.evaluate(new_mutations)

		new_cell.genome = self.genome.replicate(name=name+'g')

		# replicating other things to pass down.
		new_cell.phenotype = self.phenotype.replicate()
		new_cell.attributes = dict( self.attributes.items() )
		new_cell.parent_name = self.name
		new_cell.generation = self.generation + 1
		new_cell.cell_type = self.cell_type
		
		return new_cell

	def is_cancer(self, cancer_types = [ 2, 3 ]):
		"""
			Determines if this cell is cancerous
		"""
		return self.cell_type in cancer_types


class SelectionDistribution(object):
	@staticmethod
	def equal(cell_array):
		raw_dist = np.array([ cell.p_division() for cell in cell_array ])
		return raw_dist / float(np.sum(raw_dist))
	@staticmethod
	def cancer_only(cell_array):
		raw_dist = np.array( [ cell.p_division() if cell.is_cancer() else 0 for cell in cell_array ] )
		return raw_dist / float( np.sum(raw_dist) )


DEFAULT_CELL_ATTRIBUTES = { 'mean_mutations':150, 'chromosome_order':15, 'ploidy':2 }
class Simulator(object):
	def __init__(self, n_cells=1, cell_attributes=DEFAULT_CELL_ATTRIBUTES, phenotypes=False ):
		"""
			Creates a simulation that can be scaled
			@params:
				n_cells = 1: the number of starter cells
		"""
		cellindicies = range(n_cells)

		phenotypes = {} if not phenotypes else phenotypes

		self.cells = dict( zip( cellindicies, [ Cell( name=str(i), phenotypes=phenotypes, **cell_attributes ) for i in cellindicies ]) )
		self.attributes = {
			'cancer_created': False
		}
		self.time = 0

	def run(self, time_steps = 100, proportion_divide=0.5, stop_normal_divisions=False,**kwargs):
		"""
			Runs the model for time_steps
		"""
		auto_reduce = False

		for i in xrange(time_steps):
			if proportion_divide == 'auto_reduce' or auto_reduce:
				magnitude = kwargs.get('auto_reduce_magnitude', 1)
				auto_reduce = True
				proportion_divide = 1./float(2 + magnitude*self.time)
			num_cells_to_divide = int(proportion_divide*len(self.cells.values())+1)
			print( 'step: '+str(i+1) + ' of '+str(time_steps) + ' / Total Time: ' +str(self.time) + ' / auto_reduce: '+str(auto_reduce)+\
				', proportion_divide: '+str(proportion_divide) +' i.e. approx '+str(num_cells_to_divide)+' cells' )
			self.time += 1
			self.step( proportion_divide=proportion_divide , stop_normal_divisions=stop_normal_divisions )
	
	def create_cancer(self, update_mean_mutations=120):
		"""
			Creates a cancerous cell
		"""
		idx, celllist = zip(*self.cells.items())
		selected_cell_id = int(np.random.choice(idx, size=1, replace=False))
		selected_cancer_cell = self.cells[selected_cell_id]

		# change the cell type
		selected_cancer_cell.cell_type = 2 

		# introduce one mutation that caused it
		selected_cancer_cell.genome.mutate(number_of_new_mutations=1)
		selected_cancer_cell.genome.mean_mutations = update_mean_mutations
		selected_cancer_cell.p_division_function = pDivisionFunction.sigmoid()

		self.attributes['cancer_created'] = True
		return selected_cell_id

	def step(self, proportion_divide=0.5, stop_normal_divisions = False):

		idx, celllist = zip(*self.cells.items())

		if stop_normal_divisions:
			p_dist = SelectionDistribution.cancer_only(celllist)
		else:	
			p_dist = SelectionDistribution.equal(celllist)

		pick_size = min( int(proportion_divide*len(celllist)+1) , len(np.nonzero(p_dist)[0]) )
		cellids_to_divide = np.random.choice(idx, size=pick_size, replace=False, p = p_dist)

		biggest_index = max(idx)

		for cell in cellids_to_divide:
			biggest_index = biggest_index+1
			self.cells[biggest_index] = self.cells[cell].mitosis(name=str(biggest_index))

		# print('Took Step')
		pass