### WRITTEN BY ZAFARALI AHMED
### June 2015
## V0.2
import numpy as np

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
		"""

		self.genome_order = int ( kwargs.get( 'genome_order', 15 ) )

		self.name = kwargs.get( 'name' , '' )

		self.mutation_rate = int( kwargs.get( 'mutation_rate' , 0 ) )
		assert self.mutation_rate > -1 , ' mutation rate cannot be negative '

		self.annotations = {}
		self.mutated_loci = []
	
	def replicate ( self , name = '' ):
		"""
			replicates the current genome and all its features
			and returns the new one
			@params
				name / str
				the name of the new genome
		"""
		replicated_genome = Genome( mutation_rate = self.mutation_rate , genome_order = self.genome_order , name = name )
		replicated_genome.mutated_loci.extend( self.mutated_loci )
		replicated_genome.annotations = dict( self.annotations )

		return replicated_genome


	def mutate ( self ):
		"""
			mutates the genome and returns the loci that were mutated
		"""	
		# generate how many mutations we want
		number_of_mutations = np.random.poisson( self.mutation_rate )
		# generate random numbers representing loci


		loci = map( Mutation ,  np.around( np.random.uniform(  size = number_of_mutations ) , decimals = self.genome_order ) )
		# print loci
		# store the loci in mutated_loci if they aren't already there to represent
		
		# for locus in loci:
		# 	# if locus in self.mutated_loci:
		# 	# 	# If they are there, remove that loci to represent 
		# 	# 	self.mutated_loci.remove( locus )
		# 	# else:
		# 		# a bit flip to 1
		# 		self.mutated_loci.add( locus )

		# for locus in loci:
		# 	if not ( locus in self.mutated_loci ) :
		# 		self.mutated_loci.add( locus )
		

		## The probability of two loci being the same is extremely small 1e-10
		self.mutated_loci.extend( loci )

		# self.mutated_loci = self.mutated_loci.union( loci )

		# return mutated loci
		return loci

	def get_mutated_loci ( self , form = 'list' ):
		""" 
			@params: form / str / list
				the format of the return
			@return: list of ints
				location of the loci of the mutation (bits that are 1)
		"""
		if form == 'set':
			print 'SET is no longer supported as a method of export'
		# if form == 'set':
		# 	return self.mutated_loci

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

		# a location is given in terms of a float, create a mutation object and check if it exists
		return Mutation(location) in self.mutated_loci

	@staticmethod
	def from_mutated_loci ( mutated_loci , mutation_rate = 0 , name= '' ):
		to_return = Genome( name=name , mutation_rate = mutation_rate )
		to_return.mutated_loci = map( Mutation , sorted( list( mutated_loci ) ) ) 
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
					identification for the mutation, usually a float
				carrier / * [optional]
					the initial carrier of this mutation
					( reccomend that you store `Genome` objects )
					
		"""
		assert locus >= 0 , 'locus must be greater than or equal to 0'
		self.locus = float( locus )

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
			mutation objects are represented as #LOCI_AS_FLOAT
		"""
		return '#'+str( self.locus )

	# def __eq__ ( self, other ):
		# return self.locus == other.locus

	def __str__ ( self ):
		"""
			when converted to a string, mutation objects are just represented by their floats
		"""
		return str( self.locus )

	def to_float ( self ):
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


class GenomeCompare:
	def __init__ ( self, genomes = [ None , None ] ):
		"""
			GenomeCompare has been moved to its own file. Reimport from cc3dtools.GenomeCompare
		"""
		raise DeprecationWarning('GenomeCompare has been moved to its own file. Reimport from cc3dtools.GenomeCompare' )

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



	if method == 'naive':
		with open( file_name, 'w' ) as f:
			writer = csv.writer( f )
			for k , genome in enumerate( genomes ):
				writer.writerow( [ genome.name, genome.mutation_rate ] + sorted( genome.get_mutated_loci() ) )


