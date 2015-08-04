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
		
		g1_mutated = g1.get_mutated_loci(form='set')
		g2_mutated = g2.get_mutated_loci(form='set') 

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


		g1_mutated = g1.get_mutated_loci(form='set')
		g2_mutated = g2.get_mutated_loci(form='set')

		
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


	def genome_mutations ( self ):
		"""
			returns a list indexed by genomes
			which represent the % of genes in the 
			indexed genome which are mutated
		"""

		pass

	def draw ( self , ordering = None , save_fig = None ):

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

		if save_fig:
			plt.savefig( save_fig , format='png')
		else:
			plt.show()
		
		pass

	def get_by_name ( self , name ):
		"""
			returns a genome referenced by its @param'name'
		"""
		try:
			return self.genomes[ name ]
		except KeyError:
			try:
				return self.genomes[ str(name) ]
			except KeyError:
				try:
					return self.genomes [ int(name) ]
				except Exception as e:
					print '(!) ERROR OCCUED'+str(e)
					raise e

	@staticmethod
	def from_gen_file ( file_name , old = False ):
		"""
			imports a (unaligned) gen_file and returns a GenomeCompare object
			@params
				file_name: the name of the file containing the genomes
				old / boolean / False
					if this file was created before June 23rd 2015, it is likely to be in the float format.
					use True in that case only.
		"""

		import csv
		genomes = {}
		with open( file_name , 'r' ) as f:
			reader = csv.reader( f ) 
			for row in reader:
				if old:
					genomes[ int( row[0] ) ] = Genome.from_mutated_loci( map( float , row[2:] ) , mutation_rate = int( row[1] ) , name = int( row[0] ) ) 
				else:
					genomes[ int( row[0] ) ] = Genome.from_mutated_loci( map( int , row[2:] ) , mutation_rate = int( row[1] ) , name = int( row[0] ) ) 

		return GenomeCompare( genomes = genomes )

	@staticmethod
	def from_gen2_file( file_name , force_loci = False ):
		"""
			use a new gen_file generated from multiple chromosomes to get genomes
			@params:
				file_name: name of the file to load
				force_loci / bool / False:
					load the mutated loci into memory at run time

		"""
		from cc3dtools.Genome import gen2_to_dict
		
		genomes = gen2_to_dict( file_name , force_loci = force_loci )

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
						mapper.append( ( key , Mutation( int( mutation ) ) ) )
					mapper = dict(mapper)
					rownum += 1

				else:
					genome_name = row[1]
					mutation_rate = row[2]
					genome_mutations_rep = row[3:]
					genome_mutations_rep = map( int , map( int , genome_mutations_rep ) )
					genome_mutations = []
					for key , mutation_state in enumerate( genome_mutations_rep ):
						# print mutation_state
						if mutation_state == 1:
							genome_mutations.append( mapper[ key ] )

					genome_to_create = Genome( mutation_rate = int( mutation_rate ) , name = str( genome_name ) )
					genome_to_create.mutated_loci = genome_mutations
					genomes.append( genome_to_create )

		
		return GenomeCompare( genomes = genomes )


def newick_order( s ):
	"""
		takes in a newick representaiton of a tree
		and returns the ordering of genomes from the tree.
	"""
		# split ( brackets and then join and then split ) brakets and join. 
		# disregard the ; at the end of the string and then split according to the commas.
		# map the split strings into ints.
	return map( int , ''.join( ''.join( s.split( '(' ) ).split( ')' ) )[0:-1].split( ',' ) )
