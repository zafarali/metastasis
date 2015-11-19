import random
import time
import numpy as np
from collections import Counter

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
	
	@property
	def mutation_rate(self):
		return self.mean_mutations

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
		self.genome_order = chromosome_order
		self.ploidy = ploidy
		self.ploidy_probability = 0 # for compatibility with save_genomes2
		self.name = name
		pass
    
	@property
	def mutation_rate(self):
		return self.mean_mutations

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
		replicated_genome = Genome(mean_mutations=self.mean_mutations, chromosome_order=self.genome_order, ploidy=self.ploidy, name=name )

		replicated_chromosomes = []

		for cid, chromosome in enumerate(self.chromosomes):
			replicated_chromosomes.append( chromosome.replicate(name=name+str(cid)) )

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
			
			x = float(args[0].phenotype.get_counts()['advantageous'])
			return 1. / ( 1 + np.exp( -( a*x - c ) ) )

		return sig_fn

	@staticmethod
	def maturity(age, func=lambda: 0.5):
		"""
			Returns some division function iff the age of maturity has been reached.
		"""
		def to_return(*args, **kwargs):
			current_time = kwargs.get('time', 0)
			# dob = kwargs.get('dob', 0)
			if current_time - args[0].dob() > age:
				return func(*args, **kwargs)
			else:
				return 0

		return to_return

	@staticmethod
	def finite_divisions(number_of_divisions, func=lambda: 0.5):
		"""	
			Returns the finite division function, which returns zero
			after cell.number_of_divisions() > number_of_divisions
		"""
		def to_return(*args, **kwargs):
			print type(args[0])
			if args[0].number_of_divisions() > number_of_divisions:
				return 0
			else:
				return func(*args, **kwargs)

		return to_return
		# raise NotImplementedError('This hasn\'t been implemented yet')

# y_intercept = np.log( 1 / float( GLOBAL['_dV'] ) -1 )

# m = phenotypes[cell.id].get_counts()['advantageous']
# denominator = 1 + np.exp( -( 0.01*m - y_intercept ) )
# cell.targetVolume = min ( 1 / denominator + cell.targetVolume , GLOBAL['maxTargetVolume'] )



class Cell(object):
	def __init__(self, mean_mutations=0, chromosome_order=15, ploidy=2, name='GenericCell', cell_type=1, phenotypes=False, p_division_function=pDivisionFunction.constant(0.5), parent_name=None, max_divisions=-1):
		self.cell_type = cell_type
		self.genome = Genome(mean_mutations=mean_mutations, chromosome_order=chromosome_order, name=name+'_g', ploidy=ploidy)
		self.phenotype = Phenotype( {} if not phenotypes else phenotypes )
		self.p_division_function = p_division_function
		self.attributes = {}
		self.name = name
		self.parent_name = parent_name
		self.generation = 0
		self.date_of_birth = 0
		self.number_of_divisions = 0
		# self.ORIGINAL_DIVISON_FUNCTION = p_division_function
		self.max_divisions = max_divisions

	def p_division(self, **kwargs):
		"""
			Returns the probability of this cell dividing
		"""
		return self.p_division_function(self, **kwargs)
		# pass

	def mitosis(self, name=False, dob=0):
		"""
			Conducts a mitosis event and returns a new child cell
		"""

		# print 'Mitosis of cell_type:'+str(self.cell_type)
		if self.cell_type != 3 and self.max_divisions != -1 and self.number_of_divisions > self.max_divisions:
			# we are not a CSC
			# we do not have inifite division potential
			# and the number of divisions this cell has taken is far more than
			# the maximum number of divisions it can take.
			return None

		if not name:
			name = self.name+'child'

		self.number_of_divisions += 1

		# create new cell
		# first check if this is a CSC and :
		if self.cell_type == 3 and np.random.rand() > self.p_regeneration:
			# we are a CSC and 
			# since self.p_regeneration is small, we did not draw to create another
			# CSC cell thus we create a regular cancer cell.
			# Create a cancer cell cell that has a finite division potential.
			# with some probability it creates another one of itself.
			new_cell = Cell(name=name, cell_type = 2, p_division_function=pDivisionFunction.finite_divisions( self.max_divisions, func=self.p_division_function))
		else:
			# in all other cases, whether CSC or not
			# we just replicate ourselves
			new_cell = Cell(name=name, cell_type = self.cell_type, p_division_function=self.p_division_function)


		# pass down all traits from parent
		# in the event of CSC, we do have a max divisions
		# but since we are CSC it is bypassed in this code.
		new_cell.max_divisions = self.max_divisions

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
		new_cell.date_of_birth = dob
		
		return new_cell

	def is_cancer(self, cancer_types = [ 2, 3 ]):
		"""
			Determines if this cell is cancerous
		"""
		return self.cell_type in cancer_types

	def dob(self):
		"""	
			Returns the date of birth of this cell.
		"""
		return self.date_of_birth

	def number_of_divisions(self):
		"""
			Returns the number of divisions from birth
		"""
		return self.number_of_divisions



class SelectionDistribution(object):
	@staticmethod
	def equal(cell_array, **kwargs):
		"""
			Returns probabilities of division for each cell normalized among 
			the population
		"""
		raw_dist = np.array([ cell.p_division(**kwargs) for cell in cell_array ])
		return raw_dist / float(np.sum(raw_dist))
	@staticmethod
	def cancer_only(cell_array, **kwargs):
		"""
			Returns normalized probabilities for cancer cells only
		"""
		raw_dist = np.array( [ cell.p_division(**kwargs) if cell.is_cancer() else 0 for cell in cell_array ] )
		return raw_dist / float( np.sum(raw_dist) )
	@staticmethod
	def aged(cell_array, age=1, **kwargs):
		"""
			Returns a p_division() of choice if and only if the age is large enough.
		"""
		current_time = kwargs.get('time', 0)
		raw_dist = np.array( [ cell.p_division(**kwargs) if cell.max_divisions != 0 and cell.number_of_divisions <= cell.max_divisions else 0 for cell in cell_array ] )
		return raw_dist / float( np.sum( raw_dist ) )


DEFAULT_CELL_ATTRIBUTES = { 'mean_mutations':1, 'chromosome_order':2, 'ploidy':2 }
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
		print('Simulation started at '+str(time.ctime()))
		start_time = time.time()
		for i in xrange(time_steps):
			if proportion_divide == 'auto_reduce' or auto_reduce:
				magnitude = kwargs.get('auto_reduce_magnitude', 0.75)
				auto_reduce = True
				proportion_divide = 1./float(2 + magnitude*self.time)
			num_cells_to_divide = int(proportion_divide*len(self.cells.values())+1)
			print( 'step: '+str(i+1) + ' of '+str(time_steps) + ' / Total Time: ' +str(self.time) + ' / auto_reduce: '+str(auto_reduce)+\
				', proportion_divide: '+str(proportion_divide) +' i.e. approx '+str(num_cells_to_divide)+' cells' )
			self.time += 1
			self.step( proportion_divide=proportion_divide , stop_normal_divisions=stop_normal_divisions )
		end_time = time.time()

		print('Simulation ended at '+str(time.ctime()))
		print('Total time taken:'+str((end_time - start_time))+'seconds')

	def create_cancer(self, update_mean_mutations=2, p_division_function=pDivisionFunction.sigmoid()):
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

		# update the phenotype of the cell with a new mutation and a division function
		selected_cancer_cell.genome.mean_mutations = update_mean_mutations
		selected_cancer_cell.p_division_function = p_division_function

		self.attributes['cancer_created'] = True
		return selected_cell_id

	def create_CSC(self, update_mean_mutations=2, p_regeneration=0.02, max_divisions=4):
		cancer_id = self.create_cancer(update_mean_mutations=update_mean_mutations)
		self.cells[cancer_id].cell_type = 3
		self.cells[cancer_id].p_regeneration = p_regeneration
		self.cells[cancer_id].max_divisions = max_divisions

		return cancer_id

	def step(self, proportion_divide=0.5, stop_normal_divisions = False):

		idx, celllist = zip(*self.cells.items())

		if stop_normal_divisions:
			p_dist = SelectionDistribution.cancer_only(celllist)
		else:	
			p_dist = SelectionDistribution.equal(celllist)

		pick_size = min( int(proportion_divide*len(celllist)+1) , len(np.nonzero(p_dist)[0]) )
		cellids_to_divide = np.random.choice(idx, size=pick_size, replace=False, p = p_dist)

		# index of the biggest cell so far
		biggest_index = max(idx)
		# print 'Cells Picked: '+str(len(cellids_to_divide))+' from '+str(pick_size)
		# go over all cells that need to divide and then ask them to mitosis.
		for cell_id in cellids_to_divide:
			new_cell = self.cells[cell_id].mitosis(name=str(biggest_index), dob=self.time)
			# print 'new_cell:',new_cell
			if new_cell is not None:
				biggest_index = biggest_index+1
				self.cells[biggest_index] = new_cell


	def cell_stats(self):
		"""
			returns statistics of cells
		"""
		celllist = self.cells.values()

		cell_types = map(lambda cell: cell.cell_type, celllist)
		counts = Counter(cell_types)

		return counts

	def get_genomes(self):
		"""
			returns the genomes of cells
		"""

		return map( lambda cell: cell.genome, self.cells.values())




