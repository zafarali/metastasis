## MAY 2015

"""
	GenomeCompare
	allows easy analysis 
"""

from Genome import Genome, Mutation
import numpy as np
import matplotlib.pyplot as plt

class GenomeCompare:
	def __init__ ( self, genomes = [ None , None ] ):
		"""
			Allows comparison of multiple genomes
			@params
				genomes / dict or list
					dict must be in { genome_id : Genome() } format,
					we reccommend that genome_id must be an int
					if a list of genomes is submitted, it is converted automatically.
					to { index: Genome() }
		"""
		assert len( genomes ) > 1 , 'Must provide at least two genomes to initiatile GenomeCompare'
		# assert isinstance( genomes[0] , Genome ) and isinstance( genomes[1] , Genome ) , 'genomes must contain Genome objects only'		
		
		if type( genomes ) == list :
			# user is using an older method for genome input to maintain backward compatibility
			# we convert it to the new format and store it
			genomes = dict ( zip( range( len( genomes ) ) , genomes ) )
			print '(!) New method of inputing genomes parameter is now present, for usability, the input has been converted, refer: https://github.com/zafarali/metastasis/issues/10'

		assert type( genomes ) == dict , 'Genomes must be input in dict or list format'
		self.genomes = genomes

	def diff( self , genome1 , genome2 ):
		"""
			returns the number and loci of genes that
			are unqiue between the two genomes supplied whos id are supplied
			@params
				genome1 / int / index of the first genome to compare
				genome2 / int / index of the second genome to compare
			@return
				object with properties specified
		"""

		g1 = self.genomes[genome1]
		g2 = self.genomes[genome2]
		
		g1_mutated = set( g1.get_mutated_loci() )
		g2_mutated = set( g2.get_mutated_loci() )

		g1_only = g1_mutated.difference( g2_mutated )
		g2_only = g2_mutated.difference( g1_mutated )
		different_loci = g1_only.union( g2_only )

		different_genes = len( different_loci )

		# # check which loci in 1 are in 2
		# for locus in g1_mutated:
		# 	if not (locus in g2_mutated):
		# 		different_genes += 1
		# 		different_loci.append(locus)
		# 		g1_only.append(locus)

		# # check which loci in 2 are in 1 but not already counted
		# for locus in g2_mutated:
		# 	if not(locus in g1_mutated) and not(locus in different_loci):
		# 		different_genes +=1 
		# 		different_loci.append(locus)
		# 		g2_only.append(locus)
			
		return { 'total_unique_mutations': different_genes , 'loci_different': list( different_loci ) , 'unique_to_1' : list( g1_only ) , 'unique_to_2' : list( g2_only ) }

	def comm( self , genome1 , genome2 ):
		"""
			returns the number and loci common to
			a pair of genomes whos id are supplied
		"""

		g1 = self.genomes[genome1]
		g2 = self.genomes[genome2]


		g1_mutated = set( g1.get_mutated_loci() )
		g2_mutated = set( g2.get_mutated_loci() )

		
		common_mutation_loci = g1_mutated.intersection( g2_mutated )
		common_mutations = len( common_mutation_loci )

		# # check which loci in 1 are in 2
		# for locus in g1_mutated:
		# 	if locus in g2_mutated:
		# 		common_mutations += 1
		# 		common_mutation_loci.append(locus)

		# # check which loci in 2 are in 1 but not already counted
		# for locus in g2_mutated:
		# 	if locus in g1_mutated and not (locus in common_mutation_loci):
		# 		common_mutations += 1 
		# 		common_mutation_loci.append(locus)

		return { 'total_common_mutations': common_mutations , 'loci_common': list( common_mutation_loci ) }

	def simillarity_matrix( self ):
		"""
			returns a matrix containing the % of mutations shared
			between any two genomes
		"""

		pass

	def mutation_chart( self, draw = False, order = None ):
		"""
			returns a matrix where the y-indicies
			are genes, and 1 represents mutated
			@params
				draw / boolean / False
					if true: draws a graph representing all the
					mutations per genome

		"""
		raise DeprecationWarning('You can no longer use mutation_chart or gene_chart to draw a GC object. Use GenomeCompare.draw(...)')
		
		size = self.genomes[1].size
		num_genomes = len(self.genomes)
		rep = np.zeros( ( size, num_genomes ) )
		if order is None:
			order = [ x for x in range( num_genomes ) ]

		for genomeid in order:
			print 'genomeid:',genomeid
			genome = self.genomes[genomeid]
			for gene in genome.get_mutated_loci():
				rep[gene][genomeid] = 1

		if draw:
			plt.figure()
			plt.imshow( rep , interpolation = 'nearest' )
			plt.title('Genome Wide Visualization of Mutated Genes')
			plt.xlabel('Genome')
			plt.xticks( range( num_genomes ) , order )
			plt.ylabel('Gene')
			plt.show()
		return rep

	def gene_chart( self , draw = False ):

		return self.mutation_chart(draw=draw)

	def gene_mutations( self ):
		"""
			returns a list indexed by genes
			which represent the % of genomes over which 
			the indexed gene is mutated
		"""
		x = self.mutation_chart()
		return map ( lambda y : y / float( len( self.genomes ) ), map( sum , x ) )

	def genome_mutations ( self ):
		"""
			returns a list indexed by genomes
			which represent the % of genes in the 
			indexed genome which are mutated
		"""

		pass

	def draw ( self , ordering = None ):

		if ordering is None:
			ordering = []
			num_genomes = len(self.genomes)
			for genome_id, _ in self.genomes.items():

				ordering.append( { 
					'location':[ genome_id ], 
					'name': genome_id 
				} ) 

		plt.figure()

		for g in ordering:
			
			# get the x coordinate and plot it 
			loc = g['location'][0]
			plt.plot( [ loc ,loc ] , [ 1 , 0 ] , 'b' )

			# plot the name of the genome
			plt.text( loc , -0.002 , g['name'] , rotation = 90 , size = 'xx-small', horizontalalignment='center' )
			
			# get mutated loci and plot them
			mutated_loci = self.genomes[ int( g['name'] ) ].get_mutated_loci()

			for locus in mutated_loci: 
				plt.plot(loc, locus.to_float() , 'or')

		plt.show()
		pass

	def get_by_name ( self , name ):
		"""
			returns a genome referenced by its @param'name'
		"""

		return self.genomes[ name ]

	@staticmethod
	def from_gen_file ( file_name ):
		"""
			imports a (unaligned) gen_file and returns a GenomeCompare object
		"""

		import csv
		genomes = {}
		with open( file_name , 'r' ) as f:
			reader = csv.reader( f ) 
			for row in reader:
				genomes[ int( row[0] ) ] = Genome.from_mutated_loci( map( float , row[2:] ) , mutation_rate = int( row[1] ) , name = int( row[0] ) ) 

		return GenomeCompare( genomes = genomes )

	@staticmethod
	def from_aligned_gen_file ( file_name ):
		"""
			imports the new format of aligned gen_file and returns a GenomeCompare object
		"""
		print 'Please use GenomeCompare2 for better speeds with aligned genome files' 
		import csv
		with open ( file_name , 'r' ) as f:
			reader = csv.reader( f )
			rownum , titles = 0 , []
			mapper = []
			genomes = []
			for row in reader:
				if rownum == 0:
					mutations_in_file = row[3:]
					# create a mapper of mutation objects
					for key, mutation in enumerate( mutations_in_file ):
						mapper.append( ( key , Mutation( float( mutation ) ) ) )
					mapper = dict(mapper)
					rownum += 1

				else:
					genome_name = row[1]
					mutation_rate = row[2]
					genome_mutations_rep = row[3:]
					genome_mutations_rep = map( int , map( float , genome_mutations_rep ) )
					genome_mutations = []
					for key , mutation_state in enumerate( genome_mutations_rep ):
						# print mutation_state
						if mutation_state == 1:
							genome_mutations.append( mapper[ key ] )

					genome_to_create = Genome( mutation_rate = float( mutation_rate ) , name = str( genome_name ) )
					genome_to_create.mutated_loci = genome_mutations
					genomes.append( genome_to_create )

		
		return GenomeCompare( genomes = genomes )

		
		


class GenomeCompare2(object):
	def __init__ ( self , file_name ):
		"""
			Used to handle files that are aligned genomes. Returns a GenomeCompare2 object
			@params:
				file_name / string
					the path to the file needed to be loaded
		"""

		import csv
		with open ( file_name , 'r' ) as f:
			reader = csv.reader( f )
			rownum , titles = 0 , []
			mutation_mapper , genomes , reverse_mutation_mapper = [] , [] , [] 
			name_mapper = {}

			for row in reader:
				if rownum == 0:
					# creates a mapper
					mutations_in_file = row[3:]
					# create a mapper of key |--> mutation 
					for key, mutation in enumerate( mutations_in_file ):
						mutation_mapper.append( ( key , float( mutation ) ) ) 
						reverse_mutation_mapper.append( ( float( mutation ) , key ) )
					
					mutation_mapper = dict(mutation_mapper) # stores key --> mutation
					reverse_mutation_mapper = dict(reverse_mutation_mapper) # stores mutation --> key

					rownum += 1
				else:
					genome_name = row[1]
					
					genome_id = int(row[0])

					# create a mapper of name |--> genome_id 
					if genome_name != '':
						name_mapper[genome_name] = genome_id

					genome_mutations_rep = row[3:]
					
					genomes.append( np.array( map( lambda x: int( float( x ) ) , genome_mutations_rep ) ) ) 

		self.genomes = genomes
		self.mutation_mapper = mutation_mapper
		self.reverse_mutation_mapper = reverse_mutation_mapper
		self.name_mapper = name_mapper

	def genome_identify ( self , genome_identifier ):
		return self.name_mapper.get( genome_identifier, genome_identifier  )
	def lookup( self, genome_identifier , locus = None ):
		"""
			looks up a genome and [optional] locus
			@params
				genome_identifier / str or int
					the genome_id or name of genome you are referencing
				locus / int / None
					include if you are looking for a specific gene
		"""
		genome_id = self.genome_identify( genome_identifier )
		if locus is None:
			return self.genomes[genome_id]
		else:
			locus_idx = self.reverse_mutation_mapper.get( locus )
			return self.genomes[genome_id][locus_idx]

	def mutated_loci( self , genome_identifier ):
		genome_id = self.genome_identify( genome_identifier )
		this_genome = self.genomes[genome_id]
		non_zero = this_genome.nonzero()[0].tolist()
		to_be_returned = []
		for k in non_zero:
			to_be_returned.append( self.mutation_mapper[k] )
		return to_be_returned


	def diff( self, genome_identifier1 , genome_identifier2 ):
		id1 = self.genome_identify( genome_identifier1 )
		id2 = self.genome_identify( genome_identifier2 )

		gen1 = self.genomes[id1]
		gen2 = self.genomes[id2]

		loci = np.bitwise_xor( gen1 , gen2 )

		to_be_returned = { 'different_loci' : [] , 'total_unique_mutations' : np.sum( loci ) }

		loci = loci.nonzero()

		mutation_loci = np.vectorize( lambda locus : self.mutation_mapper[locus] )
		to_be_returned['different_loci'] = mutation_loci( loci )

		return to_be_returned

	def draw ( self , ordering ):
		plt.figure()
		for g in ordering:
			loc = g['location'][0]
			plt.plot( [ loc ,loc ] , [ 1 , 0 ] , 'b' )
			plt.text( loc , -0.002 , g['name'] , rotation = 90 , size = 'xx-small', horizontalalignment='center' )
			l = self.mutated_loci( int( g['name'] ) )
			for i in l:
				plt.plot(loc, i, 'or')
		plt.show()
		pass




def newick_order( s ):
	"""
		takes in a newick representaiton of a tree
		and returns the ordering of genomes from the tree.
	"""
		# split ( brackets and then join and then split ) brakets and join. 
		# disregard the ; at the end of the string and then split according to the commas.
		# map the split strings into ints.
	return map( int , ''.join( ''.join( s.split( '(' ) ).split( ')' ) )[0:-1].split( ',' ) )
