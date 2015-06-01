### WRITTEN BY ZAFARALI AHMED
### May 2015
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
			size / int / 1000
				number of genes/bases in genome
			genome_order / int / 4
				the order of the genome (10^genome_order)
		"""

		size = int( kwargs.get( 'size' , 1000 ) )
		assert size > 0 , 'genome_size must be non-zero positive'
		self.size = size
		self.genome_order = int ( kwargs.get( 'genome_order', 10 ) )

		self.name = kwargs.get( 'name' , '' )

		self.mutation_rate = int( kwargs.get( 'mutation_rate' , 0 ) )
		assert self.mutation_rate > -1 , ' mutation rate cannot be negative '

		self.annotations = {}
		self.mutated_loci = set()
	
	def replicate ( self ):
			replicated_genome = Genome( mutation_rate = self.mutation_rate , size = self.size , genome_order = self.genome_order )
			replicated_genome.mutated_loci = replicated_genome.mutated_loci.union( self.mutated_loci )
			replicated_genome.annotations = dict( self.annotations )

			return replicated_genome


	def mutate ( self ):
		"""self.mutated_loci = []
			mutates the genome
		"""	
		# generate how many mutations we want
		number_of_mutations = np.random.poisson( self.mutation_rate )
		# generate random numbers representing loci

		#	logic for the below compound function
		#	(1) generate uniforms
		#	(2) map that to strings
		#	(3) select the first few numbers
		#	(4) map the strings to floats
		# 	(5) map the floats to mutations
		loci = set( map( Mutation ,  np.around( np.random.uniform(  size = number_of_mutations ) , decimals = 4 ) ) )
		# print loci
		# store the loci in mutated_loci if they aren't already there to represent
		
		# for locus in loci:
		# 	# if locus in self.mutated_loci:
		# 	# 	# If they are there, remove that loci to represent 
		# 	# 	self.mutated_loci.remove( locus )
		# 	# else:
		# 		# a bit flip to 1
		# 		self.mutated_loci.add( locus )

		for locus in loci:
			if not ( locus in self.mutated_loci ) :
				self.mutated_loci.add( locus )
		
		# self.mutated_loci = self.mutated_loci.union( loci )

		# return self to allow chaining to occur
		return self

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
	def from_mutated_loci ( mutated_loci , size = 1000 , mutation_rate = 0 , genome_order = 4 ):
		to_return = Genome( size = size , genome_order = genome_order , mutation_rate = mutation_rate )
		to_return.mutated_loci = set( map( Mutation , sorted( list( mutated_loci ) ) ) )
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
		"""
		assert locus >= 0 , 'locus must be greater than or equal to 0'
		self.locus = float( locus )

	def __repr__ ( self ):
		return str( self.locus )

	# def __eq__ ( self, other ):
		# return self.locus == other.locus

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
			M = M.union( g.get_mutated_loci( form = 'set' ) )

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
				writer.writerow( sorted( genome.get_mutated_loci() ) )


