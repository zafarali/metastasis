## Author @zafarali
## June 2015
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from collections import Counter
import numpy as np
import csv
import math
import time
import pickle



def discrete_cmap(N, base_cmap='prism'):
    """ generate a shuffled discrete color map of size N based on @param base_cmap"""
 	# from: https://gist.github.com/jakevdp/91077b0cae40f8f8244a
    base = plt.cm.get_cmap(base_cmap) # get the base
    color_raw = np.linspace(0, 1, N) # generate some random colors
    np.random.shuffle(color_raw) # shuffle it
    color_list = base(color_raw) # turn them into colors
    cmap_name = base.name + str(N) # give it a random name
    return base.from_list(cmap_name, color_list, N) #return it!
 


def spatial_plot( start_file = None , end_file = None , type_colors = ( 'r', 'b', 'g' ) , format = ( 'id' , 'type' , 'x', 'y', 'z' ) , projection = '2d', hide_numbers = True ):
	"""
		displays the positions of the cellids at the time of sampling within a simulation
	"""
	assert start_file is not None and end_file is not None, 'must specify files'

	points = []

	starter_cells = [] # stores as (cellid, type)
	with open( start_file , 'r' ) as f:
		reader = csv.reader( f )
		for row in reader:
			starter_cells.append(  int( row[0] )  )


	fig = plt.figure()

	if projection == '3d':
		ax = fig.gca( projection = '3d' )

	with open( end_file , 'r' ) as f:
		reader = csv.reader( f )
		for row in reader:
			data = dict( zip( format , row ) )
			# data converstion
			data['id'] = int( data['id'] )
			data['type'] = int( data['type'] )

			data['x'] = float( data['x'] )
			data['y'] = float( data['y'] )
			data['z'] = float( data['z'] )

			data['initial'] = 1 if data['id'] in starter_cells else 0

			marker = 'x' if data['initial'] else 'o'
			if projection == '3d':
				ax.scatter( data['x'] , data['y'] , data['z'] , marker = marker , color = type_colors[ data['type'] - 1 ] )

				if not hide_numbers:
					ax.text( data['x']+0.5 , data['y']+0.5 , data['z'] , str( data['id'] ) , horizontalalignment = 'center' , color = type_colors[ data['type'] - 1 ] )
			else:
				plt.plot( data['x'] , data['y'] , marker+type_colors[ data['type'] - 1 ] )

				if not hide_numbers:
					plt.text( data['x']+0.5 , data['y']+0.5 , str( data['id'] ) , horizontalalignment = 'center' , color = type_colors[ data['type'] - 1 ] )

	if projection == '3d':
		ax.set_xlabel( 'x-axis' )
		ax.set_ylabel( 'y-axis' )
		ax.set_zlabel( 'z-axis' )
		ax.set_title('Locations of Genomes at time of final sampling', 'center')
	else:
		plt.xlabel('x-axis')
		plt.ylabel('y-axis')
		plt.title('Locations of Genomes at time of final sampling')


	plt.show()


class SpacePlot ( object ):
	def __init__ ( self , start_file = None , end_file = None , type_colors = ( 'r', 'b', 'g' ) , format = ( 'id' , 'type' , 'x', 'y', 'z' ) , projection = '2d' ):
		"""
			This class allows us to make spatial plots
		"""
		assert start_file is not None and end_file is not None, 'must specify files'

		# load up the starter cells
		starter_cells = [] # stores as (cellid, type)
		with open( start_file , 'r' ) as f:
			reader = csv.reader( f )
			for row in reader:
				starter_cells.append(  int( row[0] )  )


		cells = {}

		# Store the final positions of all the cells
		with open( end_file , 'r' ) as f:
			reader = csv.reader( f )
			for row in reader:
				data = dict( zip( format , row ) )
				# data converstion
				data['id'] = int( data['id'] )
				data['type'] = int( data['type'] )

				data['x'] = float( data['x'] )
				data['y'] = float( data['y'] )
				data['z'] = float( data['z'] )

				data['initial'] = 1 if data['id'] in starter_cells else 0

				marker = 'x' if data['initial'] else 'o'

				cells[ data['id'] ] = {
					'id': data['id'] ,
					'x': data['x'] ,
					'y': data['y'] ,
					'z': data['z'] ,
					'initial' : data['initial'] ,
					'type': data['type']
				}


		self.cells = cells
		self.starter_cells = starter_cells
		self.type_colors = type_colors
		self.format = format
		self.projection = projection
		self.start_file = start_file
		self.end_file = end_file


	def plot_all( self , hide_numbers = True ):
		"""
			plots all genomes in space according to type_colors and projection 
		"""

		# use the already built function
		# @TODO: rewrite this so we aren't re-reading the files etc.
		spatial_plot( start_file = self.start_file , \
			end_file = self.end_file , \
			type_colors = self.type_colors , \
			format = self.format , \
			projection = self.projection , \
			hide_numbers = hide_numbers )
		pass


	def plot_selected( self , *args , **kwargs ):
		"""
			plots the cells in clusters whos ids are contained in lists
			within *args according to a random color map
			eg: args = ( [1,2,3,4], [5,6,7,8] )
			this will print 1,2,3,4 in one color and 5,6,7,8 in another
			@params:
				args:	[mandatory]
					sequence of clusters of Individuals 
				depth / str / 'UNKNOWN DEPTH': 
					depth of the plot
				hide_numbers / bool / True
					hide the display of numbers in the plot


				*these parameters can be overriden for this plot:
				type_colors: 
					color mappings for each type of individual 
				projection:
					2d or 3d

		"""

		projection = kwargs.get('projection', self.projection)
		type_colors = kwargs.get('type_colors', self.type_colors)
		title = ' Locations of Genomes at time of final sampling lineage depth: ' + str( kwargs.get('depth', 'UNKNOWN') )
		hide_numbers = kwargs.get( 'hide_numbers' , True )
		fig = plt.figure()

		if projection == '3d':
			ax = fig.gca( projection = '3d' )

		colormap = discrete_cmap( 900 )
		num_args = float ( len( args ) )

		plotted_ids = []

		# this function converts an individual to cell_data
		# it also accounts for the fact that some cells might have
		# died and not made it into the final sampling
		def individual_to_cell_data( individual ):
			try:
				return self.cells[ individual.id ]
			except KeyError:
				return {
					'x':0,
					'y':0,
					'z':0,
					'initial':-1, # -1 represents that the cell died.
					'id': individual.id,
				}

		# first plot the clusters
		for cluster_id, cluster in enumerate( args ):
			# we map the individuals id in the cluster to obtain all the data
			cells_with_data = map( individual_to_cell_data , cluster )
			
			for data in cells_with_data:

				# cell died, thus we just skip plotting it
				if data['initial'] == -1:
					continue

				# add the id of the cell to the plotted_ids so we plot everything else in grey
				plotted_ids.append( data['id'] )

				selected_color = colormap( cluster_id )
				marker = 'x' if data['initial'] else 'o'

				if projection == '3d':
					ax.scatter( data['x'] , data['y'] , data['z'] , marker = marker , color = selected_color )

					if not hide_numbers:
						ax.text( data['x']+0.5 , data['y']+0.5 , data['z'] , str( data['id'] ) , horizontalalignment = 'center' , color = selected_color )
				else:
					plt.plot( data['x'] , data['y'] , marker = marker, color = selected_color )
					
					if not hide_numbers:
						plt.text( data['x']+0.5 , data['y']+0.5 , str( data['id'] ) , horizontalalignment = 'center' , color = selected_color )
		



		# everything else must be plotted in gray!
		for _ , data in self.cells.items():
			
			# skip if this has already been plotted!
			if data['id'] in plotted_ids: 
				continue

			# select the color gray
			selected_color = (.5,.5,.5,1.0)


			marker = 'x' if data['initial'] else 'o'
			if projection == '3d':
				ax.scatter( data['x'] , data['y'] , data['z'] , marker = marker , color = selected_color )

				if not hide_numbers:
					ax.text( data['x']+0.5 , data['y']+0.5 , data['z'] , str( data['id'] ) , horizontalalignment = 'center' , color = selected_color )
			else:
				plt.plot( data['x'] , data['y'] , marker = marker, color = selected_color )

				if not hide_numbers:
					plt.text( data['x']+0.5 , data['y']+0.5 , str( data['id'] ) , horizontalalignment = 'center' , color = selected_color )

		if projection == '3d':
			ax.set_xlabel( 'x-axis' )
			ax.set_ylabel( 'y-axis' )
			ax.set_zlabel( 'z-axis' )
			ax.set_title(title, 'center')
		else:
			plt.xlabel('x-axis')
			plt.ylabel('y-axis')
			plt.title(title)



		plt.show()

		pass

class PostProcess( object ):
	def __init__ ( self , end_file , gc , format = ( 'id' , 'type' , 'x', 'y', 'z' ) , pickle_import = False  ):
		"""
			allows us to get pairwise distances, shared and private mutations
			between all cells.
			@params
				end_file / str / [mandatory]
				gc / GenomeCompare / [mandatory]
		"""
		if not pickle_import:
			self.gc = gc 
			self.cell_locations = {}
			self.cells_available = []
			# save the locations of all the cells
			# also save their types
			with open( end_file, 'r' ) as f:
				reader = csv.reader( f )

				for row in reader:

					data = dict ( zip ( format , row ) )

					data['id'] = int( data['id'] )
					data['type'] = int( data['type'] )

					data['x'] = float( data['x'] )
					data['y'] = float( data['y'] )
					data['z'] = float( data['z'] )

					self.cells_available.append( data['id'] ) # save available cell ids
					self.cell_locations[ data['id'] ] = ( data['x'] , data['y'] , data['z'] , data['type'] ) # save tuple of details

			self.__frompickle__ = False
			self.__executed__ = False


	def execute( self ):
		"""
			Calculates the number of shared and private mutations and the distance between every cell
			(!) warning, this is a computationally heavy function

		"""
		assert self.__executed__ == False, 'Cannot re-execute an already executed PostProcess object'
		pass

		results = {}
		
		num_loops = 0

		start_time = time.time()

		# tuple information
		dist_vs_shared = []
		dist_vs_private = []
		dist_vs_proportion = []

		# individual lists with relevant information in corresponding indicies

		_distances = []
		_shared = []
		_private = []
		_type2 = []
		_type1 = []

		print 'start time:',start_time
		for index, i in enumerate(self.cells_available):

			results[i] = {}

			for j in self.cells_available:
				
				num_loops += 1

				if i == j: # do not compute when i == j as they are just 0
					shared = len( self.gc.genomes[i].get_mutated_loci() )
					private = 0
					distance = 0
					prop = shared / float( private + shared )


					results[i][j] = { 'distance' : distance , 'private' : private , 'shared' : shared }
					results[j][i] = { 'distance' : distance , 'private' : private , 'shared' : shared }


					dist_vs_private.append( ( distance , private , self.cell_locations[i][3] , self.cell_locations[j][3] ) )
					dist_vs_shared.append( ( distance , shared , self.cell_locations[i][3] , self.cell_locations[j][3] ) )
					dist_vs_proportion.append( ( distance , prop , self.cell_locations[i][3] , self.cell_locations[j][3] ) )

					# specific lists
					_distances.append( distance )
					_shared.append( shared )
					_private.append( private )
					_type1.append( self.cell_locations[i][3] )
					_distances.append( self.cell_locations[j][3] )


					# we reach the middle of the matrix, thus we break and skip to the next coloumn.
					# print 'completed ',i
					break

				a = self.cell_locations[i]
				b = self.cell_locations[j]

				distance = distance_between( a , b )
				shared = self.gc.comm( i , j )['total_common_mutations']
				private = self.gc.diff( i , j )['total_unique_mutations']

				prop = shared / float( private + shared )


				# store distance vs. private and shared tuples for easy plotting later
				dist_vs_private.append( ( distance , private , a[3] , b[3] ) )
				dist_vs_shared.append( ( distance , shared , a[3] , b[3] ) )
				dist_vs_proportion.append( ( distance , prop , a[3] , b[3] ) )

				_distances.append( distance )
				_shared.append( shared )
				_private.append( private )
				_type1.append( a[3] )
				_type2.append( b[3] )

				# the resulting matrix is going to be symmetric
				results[i][j] = { 'distance' : distance , 'private' : private , 'shared' : shared }
				results[j][i] = { 'distance' : distance , 'private' : private , 'shared' : shared }

			# end for j
		# end for i


		self._distances = np.array( _distances )
		self._shared = np.array( _shared )
		self._private = np.array( _private )
		self._type1 = np.array( _type1 )
		self._type2 = np.array( _type2 )
		self._proportion = self._shared / np.float_( self._private + self._shared )

		self.data = results

		self.dist_vs_shared = dist_vs_shared
		self.dist_vs_private = dist_vs_private
		self.dist_vs_proportion = dist_vs_proportion

		self.__executed__ = True

		end_time = time.time()
		
		print 'Completed\n Total Number of comparisons:' + str( num_loops )
		print 'Number of cells ' + str( len( self.cells_available ) )
		print 'Total Time:' + str( end_time - start_time ) + 's'

	def pwise ( self , c1 , c2 ):
		"""
			Returns pairwise information for cells c1 and c2
		"""
		assert self.__executed__ == True, 'You must first PostProcess.execute() before you can access other methods'

		return self.data[c1][c2]

		pass

	def frequency_analyze( self , cellids ):
		"""
			Returns a tuple with the frequency spectra of the ( # of shared mutations ) VS ( # of cells in cluster that share it )
			(!) warning, this is a computationally heavy function
			@params:
				cellids / list of int / [mandatory]
					the list of cell ids that you want to sample from the larger space
					for best results and to actually get 'cluster' data you will select cellids that are close together 
					spatially.
		"""
		
		counter = Counter()

		for cellid in cellids:
			# get the loci 
			mutated_loci = self.gc.get_by_name( cellid ).get_mutated_loci( form = 'set' )

			# counts the number of times a given mutation appears in the cluster
			counter.update( mutated_loci )

		# counts the number of a count appears in the cluster
		return Counter( [ v for _, v in counter.most_common() ] )

	def tumor_COM( self ):
		"""
			calculates the COM of the tumor
		"""

		# get all cancer cells, extract the location data from it into a numpy array
		r_vectors = map( lambda x: np.array([ x[1][0], x[1][1], x[1][2] ]) , filter( lambda x: x[1][3] == 2 or x[1][3] == 3, self.cell_locations.items() ) )

		# since it is now in a numpy array, we can do elementwise-summation:
		R_CM = np.sum( r_vectors , axis = 0 ) / float( len( r_vectors ) )

		return { 'x' : R_CM[0] , 'y' : R_CM[1] , 'z' : R_CM[2] }


	def pickle_save( self ):
		raise FutureWarning('This function is yet to be implemented')
		assert self.__executed__ == True, 'You must first PostProcess.execute() before you can access other methods'

		pickle.dump( self.data , file( 'pp.executed.pickle' , 'w' ) )

	@staticmethod
	def from_pickle( file_name = 'pp.executed.pickle' ):
		raise FutureWarning('This function is yet to be implemented')
		to_return = PostProcess( 'x' , 'x' , pickle_import = True )

		to_return.data = pickle.load( file( file_name , 'r' ) )
		to_return.__executed__ = True
		to_return.__frompickle__ = True



def distance_between ( a , b ):
	q = b[0] - a[0] 
	p = b[1] - a[1] 
	r = b[2] - a[2]
	return math.sqrt( (p * p) + (q * q) + (r * r) )




