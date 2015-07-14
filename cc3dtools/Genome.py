### WRITTEN BY ZAFARALI AHMED
### July 2015
## V0.3
import numpy as np

class Chromosome(object):
	def __init__ ( self , **kwargs ):
		"""
			holds a chromosome
			@params:
			mutation_rate / int / 0
				the mean number of new mutations that will be produced
				every mutate() event.
			chromosome_order / int / 15
				the number of loci in the genome ( 10^chromosome_order )
			name / str / ''
				the name of the chromosome.
		"""

		self.name = kwargs.get( 'name' , '' )

		self.chromosome_order = int ( kwargs.get( 'chromosome_order' , 15 ) )
		self.mutation_rate = int ( kwargs.get( 'mutation_rate' , 0 ) )
		assert self.mutation_rate > - 1 , 'mutation rate cannot be negative'

		self.annotations = {}
		self.mutated_loci = set()

	def replicate( self , name = '' ):
		"""
			replicates the current chromosome and all its features
			and returns a new one
			@params
				name / str
					the name of the new chromosome
		"""

		replicated_chromosome = Chromosome( mutation_rate = self.mutation_rate , chromosome_order = self.chromosome_order , name = name )
		replicated_chromosome.mutated_loci = replicated_chromosome.mutated_loci.union( self.mutated_loci )

		replicated_chromosome.annotations = dict( self.annotations )

		return replicated_chromosome

	def mutate( self ):
		"""
			Mutates the chromosome
			@returns:
				set / the loci which were mutated
		"""

		# generate the number of mutations needed
		number_of_mutations = np.random.poisson( self.mutation_rate )
		number_of_loci = 10**self.chromosome_order

		loci = set( map( Mutation , np.random.randint( number_of_loci , size = number_of_mutations ) ) )


		self.mutated_loci = self.mutated_loci.union( loci )

		# return the newly mutated loci / can be used by cc3dtools.Phenotype
		return loci

	def get_mutated_loci ( self , form = 'list' ):
		""" 
			@params: form / str / list
				the format of the return
			@return: list of ints
				location of the loci of the mutation (bits that are 1)
		"""

		if form == 'set':
			return self.mutated_loci

		return list( self.mutated_loci )

	def annotate( self , locus , name ):
		self.annotations[name] = locus

	def is_mutated ( self , **kwargs ):
		"""
			is_mutated
				returns if a locus or specified annotation is mutated (i.e bit == 1)
			@params
				locus: index of gene
				name: annotation of the locus (use Genome.annotate to annotate loci)
				annotation: annotation of the locus (use Genome.annotate to annotate loci)
			@return
				boolean: depending on if the gene has been mutated
		"""

		
		# get locus
			# if no locus, get the name
			# get the locus of that name
			# if no locus related to the name, return None

		location = kwargs.get( 'locus' , self.annotations.get( kwargs.get( 'name' , None ) or kwargs.get( 'annotation' , None ) , None ) )

		# if no location is given return false
		if location is None:
			return False

		# if a mutation is given as a search term, check for the mutation itself
		if isinstance( location , Mutation ):
			return location in self.mutated_loci

		# a location is given in terms of a int, create a mutation object and check if it exists
		return Mutation(location) in self.mutated_loci
	
	@staticmethod
	def from_mutated_loci ( mutated_loci , mutation_rate = 0 , name = '' ):
		to_return = Chromosome( name = name , mutation_rate = mutation_rate )
		to_return.mutated_loci = set( map( Mutation , sorted( list( mutated_loci ) ) ) )

		return to_return



"""
	Genome
		emulates a 'genome' of a cell. has replication and mutation capabilities
		the genome contains 0s and 1s. If a base is mutated it is represented by 1
"""
class Genome(object):
	def __init__ ( self, **kwargs ):
		"""
			creates a new genome
			@params:
			mutation_rate / float or int / 0
				the rate at which bases of the genome are mutated (i.e bit flipped from 0 --> 1)
			genome_order / int / 10
				the order of the genome (10^genome_order)
			name / str / ''
				name of the genome
			ploidy / int / 2
				number of sets of chromosomes
			ploidy_probability / float / 0
				the probability of a ploidy event
		"""

		self.genome_order = int ( kwargs.get( 'genome_order', 15 ) )

		self.name = kwargs.get( 'name' , '' ) 
		
		self.ploidy_probability = float ( kwargs.get( 'ploidy_probability' , 0 ) )
		assert self.ploidy_probability <= 1 and self.ploidy_probability >= 0, 'ploidy_probability must be between 0 and 1'

		ploidy = kwargs.get( 'ploidy' , 1 )

		self.mutation_rate = int( kwargs.get( 'mutation_rate' , 0 ) )
		assert self.mutation_rate > -1 , ' mutation rate cannot be negative '

		self.annotations = {}
		

		self.chromosomes = [ Chromosome( mutation_rate = self.mutation_rate , chromosome_order = self.genome_order , name = self.name + '_' + str(i) ) for i in xrange( ploidy ) ]

		self.mutated_loci = list()

	
	def replicate ( self , name = '' ):
		"""
			replicates the current genome and all its features
			and returns the new one
			@params
				name / str
				the name of the new genome
		"""

		replicated_chromosomes = []


		for i , chromosome in enumerate( self.chromosomes ):
			replicated_chromosomes.append( chromosome.replicate( name = chromosome.name ) )

			replication_error = np.random.uniform() < self.ploidy_probability

			if replication_error:
				if np.random.uniform() < 0.5:
					# replicate this chromosome twice
					replicated_chromosomes.append( chromosome.replicate( name = chromosome.name+'_replicated' ) )
				else:
					# delete a random chromsome
					replicated_chromosomes.pop( np.random.randint( len(replicated_chromosomes) ) )


		replicated_genome = Genome( mutation_rate = self.mutation_rate , genome_order = self.genome_order , name = name )
		
		replicated_genome.chromosomes = replicated_chromosomes

		return replicated_genome


	def mutate ( self , unique_only = False):
		"""
			mutates the genome and returns the loci that were mutated
		"""	

		loci = []

		for chromosome in self.chromosomes:
			loci.extend( chromosome.mutate() )		

		# update the mutated_loci array
		self.mutated_loci.extend( loci )

		if unique_only:
			return set(loci)

		# return mutated loci
		return loci


	def get_mutated_loci ( self , form = 'list', method = 'list', unique_only = False ):
		""" 
			@params: form / str / list
						the format of the return
					method / str / list
						the format of the return
					unique_only / bool / False
						returns only unique mutations as a set()

			@return: list of ints
					location of the loci of the mutation (bits that are 1)
		"""
		# if form == 'set':
		# 	print 'SET is no longer supported as a method of export'

		if form == 'set' or method == 'set' or unique_only:
			return set(self.mutated_loci)

		return self.mutated_loci 

	def annotate ( self , locus , name ):
		self.annotations[name] = locus

	def is_mutated ( self , **kwargs ):
		"""
			is_mutated
				returns if a locus or specified annotation is mutated (i.e bit == 1)
			@params
				locus: index of gene
				name: annotation of the locus (use Genome.annotate to annotate loci)
				annotation: annotation of the locus (use Genome.annotate to annotate loci)
			@return
				boolean: depending on if the gene has been mutated
		"""
		raise DeprecationWarning('Genome.is_mutated is no longer supported')

	@staticmethod
	def from_chromosome_data ( chromosome_data , mutation_rate = 0 , name = '' , ploidy_probability = 0 , genome_order = 15 ):
		"""
			Generates a Genome from an array of chromsome loci
		"""
		to_return = Genome( mutation_rate = mutated_loci , ploidy_probability = ploidy_probability, genome_order = genome_order )
		chromosomes = []
		mutated_loci = []
		for chromosome in chromosome_data:
			chromsomes.append( Chromosome.from_mutated_loci( chromosome['loci'] , mutation_rate = chromosome['mutation_rate'], name= chromosome['name']) )
			mutated_loci.extend( chromsome['loci'] )

		to_return.mutated_loci = mutated_loci
		to_return.chromosomes = chromosomes

		return to_return


"""
	Mutation
		stores a mutation at a single locus
"""
class Mutation(object):
	def __init__ ( self , locus , **kwargs ):
		"""
			@params:
				locus / *
					identification for the mutation, usually a int
				carrier / * [optional]
					the initial carrier of this mutation
					( reccomend that you store `Genome` objects )
					
		"""
		assert locus >= 0 , 'locus must be greater than or equal to 0'

		# for backwards compatibility
		if type( locus ) is float:
			locus = locus * ( 10 ** 15 )

		self.locus = int( locus )

		initial_carrier = kwargs.get( 'carrier' , None )
		self.carriers = [ initial_carrier ] if initial_carrier else []

	def add_carrier( self , carrier ):
		"""
			will add a carrier to this object
			@params:
				carrier / * 
					carrier of this mutation to add
					( reccomend that you store `Genome` objects )
		"""
		self.carriers.append( carrier )
		pass

	def __repr__ ( self ):
		"""
			mutation objects are represented as #LOCI
		"""
		return '#'+str( self.locus )

	# def __eq__ ( self, other ):
		# return self.locus == other.locus

	def __str__ ( self ):
		"""
			when converted to a string, mutation objects are just represented by their loci
		"""
		return str( self.locus )

	def to_float ( self ):
		num_digits = len( str( self.locus ) )
		return self.locus / 10**float( num_digits )

	def to_int ( self ):
		return self.locus


	def __cmp__ ( self , other ):
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


def save_genomes2( genomes , file_name = 'genomes_saved_output'):
	"""
		method for saving new genomes, supply an array of Genome objects
	"""
	with open( file_name + '.gen2' , 'w' ) as f:
		writer = csv.writer( f )
		for genome in genomes:
			writer.writerow( [ 'G', genome.name, genome.mutation_rate, genome.ploidy_probability , genome.genome_order ] )
			for chromsome in genome.chromsomes:
				writer.writerow( [ chromosome.name, chromsome.mutation_rate ] + sorted( chromsome.get_mutated_loci() ) )
	print('saved genome data to ' + file_name + '.gen2' )

def save_genomes( genomes , file_name = 'genomes_saved_output.csv' , method = 'naive' ):
	"""
		saves an array of genomes to a file
		@params:
			genomes / array
				an array of genomes to be genomes to be saved
			file_name / string / 'genomes_saved_output.csv'
				file name to save the genomes into
			method / string / 'naive'
			 'naive'
				the method of saving the genome (naive saves it in the old fashion / rows correspond to genomes)
			 'aligned'
			 	stores the loci in the new advanced format

		(!) this only saves the mutated loci and not genome sizes etc.
	"""
	assert genomes is not None ,  'you must supply an array of genomes as the first argument'
	assert len( genomes ) > 0 , 'you must supply at least one genome'
	import csv
	
	if method == 'aligned':

		M = set() # all unique mutated loci will be stored here
		for g in genomes:
			M = M.union( set( g.get_mutated_loci( form = 'set' ) ) )

		M = sorted( list( M ) ) # this is the order in which mutations will be saved

		storage = np.zeros( ( len( genomes ) , len( M ) ) )

		# create a mapping
		mapping = []
		for k , m in enumerate( M ):
			mapping.append( ( m , k ) )

		mapping = dict( mapping )

		for gid , g in enumerate( genomes ):
			mutations = sorted( g.get_mutated_loci() )
			
			for locus in mutations:
				indx = mapping[locus]
				# print 'mapping '+str(locus)+'to '+str(indx)
				storage[ gid , indx ] = 1

		titles = [ 'genomeid', 'name', 'mutation_rate' ] + M


		with open( file_name , 'w' ) as f:
			writer = csv.writer( f )
			writer.writerow( titles )
			for rowid , row in enumerate( storage ):
				writer.writerow( [ rowid , genomes[rowid].name , genomes[rowid].mutation_rate ] + list( row ) )
	print 'saved genome data to '+file_name


	if method == 'naive':
		with open( file_name, 'w' ) as f:
			writer = csv.writer( f )
			for k , genome in enumerate( genomes ):
				writer.writerow( [ genome.name, genome.mutation_rate ] + sorted( genome.get_mutated_loci() ) )
		print 'saved genome data to '+file_name

def csv_to_gen2( file_name ):
	"""
		converts an old genome saved in the csv format into the new gen2 format
	"""
	pass
