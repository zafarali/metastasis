## Author @zafarali
## June 2015
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.path import Path
from mpl_toolkits.mplot3d import Axes3D
from collections import Counter
import numpy as np
import csv
import math
import time
import pickle



def random_color_map(N, base_cmap='prism'):
    """ generate a shuffled discrete color map of size N based on @param base_cmap"""
 	# from: https://gist.github.com/jakevdp/91077b0cae40f8f8244a
    base = plt.cm.get_cmap(base_cmap) # get the base
    color_raw = np.linspace(0, 1, N) # generate some random colors
    np.random.shuffle(color_raw) # shuffle it
    color_list = base(color_raw) # turn them into colors
    cmap_name = base.name + str(N) # give it a random name
    return base.from_list(cmap_name, color_list, N) #return it!

def discrete_cmap(N, base_cmap='prism'):
	""" generate a discrete color map of size N based on @param base_cmap"""
	# from: https://gist.github.com/jakevdp/91077b0cae40f8f8244a
	if N == 1:
		return plt.cm.get_cmap('jet')
	base = plt.cm.get_cmap(base_cmap)
	color_list = base(np.linspace(0,1,N+1))
	# color_list[0] = (.5,.5,.5,1.0)
	# print color_list
	cmap_name = base.name+str(N)
	return base.from_list(cmap_name, color_list, N+1)


def spatial_plot( start_file = None , end_file = None , type_colors = ( 'r', 'b', 'g' ) , format = ( 'id' , 'type' , 'x', 'y', 'z' ) , projection = '2d', hide_numbers = True , plot_stack = None , save_fig = None):
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

	plt.axes().set_aspect('equal', 'datalim')
	if plot_stack:
		while len( plot_stack ):
			args = plot_stack.pop()
			plt.plot(*args)

	if save_fig:
		plt.savefig( save_fig , format='png')
		plt.clf()
		plt.cla()
	else:
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


	def plot_all( self , hide_numbers = True , plot_stack = None , save_fig = None ):
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
			hide_numbers = hide_numbers , \
			plot_stack = plot_stack , \
			save_fig = save_fig )
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
		save_fig = kwargs.get( 'save_fig', None )
		type_colors = kwargs.get('type_colors', self.type_colors)
		title = ' Locations of Genomes at time of final sampling lineage depth: ' + str( kwargs.get('depth', 'UNKNOWN') )
		hide_numbers = kwargs.get( 'hide_numbers' , True )
		fig = plt.figure()

		if projection == '3d':
			ax = fig.gca( projection = '3d' )

		colormap = random_color_map( len(args) )
		num_args = float ( len( args ) )

		plotted_ids = []

		# this function converts an individual to cell_data
		# it also accounts for the fact that some cells might have
		# died and not made it into the final sampling
		def individual_to_cell_data( individual ):
			if type(individual) is int:
				return self.cells[ individual ]
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
			except AttributeError:
				return {
					'x':0,
					'y':0,
					'z':0,
					'initial':-1,
					'id':individual
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


		plt.axes().set_aspect('equal', 'datalim')

		if save_fig:
			plt.savefig( save_fig , format='png')
			plt.clf()
			plt.cla()
		else:
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

	def frequency_analyze( self , cellids , return_loci = False ):
		"""
			Returns a tuple with the frequency spectra of the ( # of shared mutations ) VS ( # of cells in cluster that share it )
			@params:
				cellids / list of int / [mandatory]
					the list of cell ids that you want to sample from the larger space
					for best results and to actually get 'cluster' data you will select cellids that are close together 
					spatially.
			@returns:
				Frequency / Counter and the number of cells / int
		"""
		
		counter = Counter()

		for cellid in cellids:
			# get the loci 
			try:
				mutated_loci = self.gc.get_by_name( cellid ).get_mutated_loci( form = 'set' )
			except Exception as e:
				print cellid, cellids
				print str(e)
				raise e


			# counts the number of times a given mutation appears in the cluster
			counter.update( mutated_loci )

		if return_loci:
			return counter

		# counts the number of a count appears in the cluster, return the total number of cells in that cluster
		return Counter( [ v for _, v in counter.most_common() ] ) , len( cellids )

	def tumor_COM( self ):
		"""
			calculates the COM of the tumor
			@returns:
				tuple containing the COM
				( x , y , z )
		"""

		# get all cancer cells, extract the location data from it into a numpy array
		r_vectors = map( lambda x: np.array([ x[1][0], x[1][1], x[1][2] ]) , filter( lambda x: x[1][3] == 2 or x[1][3] == 3, self.cell_locations.items() ) )

		# since it is now in a numpy array, we can do elementwise-summation:
		R_CM = np.sum( r_vectors , axis = 0 ) / float( len( r_vectors ) )

		return ( R_CM[0] , R_CM[1] , R_CM[2] )

	def nearest( self , x , y , z , radius = 5 , type_restrictions = None ):
		"""
			returns the cellids and the x,y,z coordinates of the cells within radius
			of x, y, z
			@params:
				x,y,z / float,float,float
					coordinates around which we look for nearest cells
				radius / int / 5
					radius around which we must search
				type_restrictions / list / None
					the types to which we restrict our sampling
			@return:
				list of dicts of the results of the form:
					[ { 'id': id, 'x': x, 'y':y, 'z':z, 'type':type } , ... ]
		"""

		filtered_list = self.cell_locations.items()

		if type_restrictions:
			assert type( type_restrictions ) is list , 'type_restrictions must be a list of ints representing types'
			filtered_list = filter( lambda x: x[1][3] in type_restrictions, filtered_list )

		# first map the cell_locations into a new dict
		r_vectors = map( lambda x: { 'id': x[0], 'x': x[1][0], 'y': x[1][1], 'z': x[1][2], 'type':x[1][3] } ,  filtered_list )
		return filter( lambda r: ( r['x'] - x )**2 + ( r['y'] - y )**2 + ( r['z'] - z )**2 <= radius**2, r_vectors )

	def cells_in_ellipse_at( self , x , y , z , radii , type_restrictions = None ):
		filtered_list = self.cell_locations.items()

		if type_restrictions:
			assert type( type_restrictions ) is list , 'type_restrictions must be a list of ints representing types'
			filtered_list = filter( lambda x: x[1][3] in type_restrictions, filtered_list )

		# first map the cell_locations into a new dict
		r_vectors = map( lambda x: { 'id': x[0], 'x': x[1][0], 'y': x[1][1], 'z': x[1][2], 'type':x[1][3] } ,  filtered_list )
		return filter( lambda r: ( ( r['x'] - x ) / radii[0] )**2 + ( (  r['y'] - y  )/radii[1] )**2 + ( ( r['z'] - z )/radii[2])**2 <= 1, r_vectors )

	def cluster_return( self , *args, **kwargs ):
		"""	
			DEPRECATED use PostProcess.sample_circular()
		"""
		raise DeprecationWarning('PLEASE USE PostProcess.sample_circular() in future code.')

	def sample_circular( self , x , y , z , theta, step_size , steps , cluster_radius , **kwargs ):
		"""
			searches in incremental step_size's from x,y,z and returns the nearest neighbours in the circular 
			regions of cluster_radius, we travel in a theta direction.
			@params:
				x,y,z / float,float,float 
					location from which we want to start our search
				theta / int
					angle (in radians) at which we want to search
				step_size / int
					the step sizes we want to increment our search by
				steps / int 
					the total number of steps to take
				cluster_radius / int
					the radius of the cluster we wish to sample
				type_restrictions / list / None
					restrict the sampling to cells of certain types
				show_line_plot / bool / False
					draw the sampling on a graph
				return_plot_stack / bool
					returns a plot stack for drawing

		"""
		return self.sample_ellipsoid( x, y, z, theta, step_size, steps, ( cluster_radius, cluster_radius, 1 ), **kwargs )

	def sample_ellipsoid( self , x , y , z , theta , step_size , steps , radii , type_restrictions = None , show_line_plot = False , return_plot_stack = False ):
		"""
			searches in incremental step_size's from x,y,z and returns the nearest neighbours in the circular 
			regions of cluster_radius, we travel in a theta direction.
			@params:
				x,y,z / float,float,float 
					location from which we want to start our search
				theta / int
					angle (in radians) at which we want to search
				step_size / int
					the step sizes we want to increment our search by
				steps / int 
					the total number of steps to take
				radii / tuple / ( int , int )
					a tuple of (semiminor_axis , semimajor_axis) for scaling the ellipse
				type_restrictions / list / None
					restrict the sampling to cells of certain types
				show_line_plot / bool / False
					draw the sampling on a graph

		"""

		sin_theta = np.sin( theta )
		cos_theta = np.cos( theta )

		results = []
		plot_stack = []
		# create the circle template
		if show_line_plot or return_plot_stack:
			pnts = np.linspace( 0 , 2 * np.pi , 100 )
			circle_x = radii[0] * np.sin( pnts )
			circle_y = radii[1] * np.cos( pnts )


		for step in range( steps ):

			# calculate the position of the sampling
			distance_travelled = step * step_size
			position_x = x + distance_travelled * cos_theta
			position_y = y + distance_travelled * sin_theta

			# do the sampling
			sample = self.cells_in_ellipse_at( position_x , position_y , z , radii = radii , type_restrictions = type_restrictions )
			sample_cellids = [ cell['id'] for cell in sample ]

			#analyze that sample
			results.append( ( distance_travelled , sample_cellids ) )

			if show_line_plot:
				plt.plot( x + circle_x + distance_travelled * cos_theta , y + circle_y + distance_travelled * sin_theta)

			if return_plot_stack:
				plot_stack.append( [ x + circle_x + distance_travelled * cos_theta , y + circle_y + distance_travelled * sin_theta ] )
		#endfor

		if show_line_plot:
			plt.plot( [x, x + step_size * steps * cos_theta], [ y, y + step_size * steps * sin_theta] )

		if return_plot_stack:
			plot_stack.append( [ [ x, x + step_size * steps * cos_theta], [ y, y + step_size * steps * sin_theta] ]  )

		if return_plot_stack:
			return results, plot_stack
		else:
			return results

		pass

	def sample_rectangle ( self , x , y , z , theta , step_size , steps , edge_lengths, **kwargs ):
		assert len(edge_lengths) > 1, 'need to have more than two edges'

		rect_vertices = [ [ 0 , 0 ] , [ 0 , edge_lengths[0] ], [ edge_lengths[1] , edge_lengths[0] ] , [ edge_lengths[1] , 0 ] ]

		return self.sample_polygon( x, y, z, theta, step_size, steps, rect_vertices , **kwargs )
	
	def sample_square ( self , x , y , z , theta , step_size , steps , width , **kwargs ):

		square_vertices = [ [ 0 , 0 ] , [ 0 , width ], [ width , width ] , [ width , 0 ] ]

		return self.sample_polygon( x, y, z, theta, step_size, steps, square_vertices , **kwargs )

	def sample_polygon ( self , x , y , z , theta , step_size , steps , polygon_vertices , type_restrictions = None , show_line_plot = False , return_plot_stack = False ):
		"""
			samples polgons defined by the polygon_vertices, starting from (x,y,z) along the theta direction
			for steps of step_size
			@params:
				x,y,z / ints or floats
					location of the start of the sampling procedure
				theta / ints or floats
					direction along which sampling should occur
				step_size / int
					the steps after which sampling should be conducted
				steps / int
					the number of samples to take
				polygon_vertices / tuple of ints or floats
					the coordinates of the vertices of the polygon (should be relative to the origin)
				type_restrictions / tuple of ints / None
					the types to which sampling should be restricted
				show_line_plot / bool / False
					plot the sketch of the path
				return_plot_stack / bool / False
					the plot stack will be returned for future plotting in cc3dtools.SpacePlot
		"""

		filtered_list = self.cell_locations.items()
		
		# remove the points which do not have the required type.
		if type_restrictions:
			assert type( type_restrictions ) is list , 'type_restrictions must be a list of ints representing types'
			filtered_list = filter( lambda x: x[1][3] in type_restrictions, filtered_list )

		# obtain cell information
		r_vectors = map( lambda x: { 'id': x[0], 'x': x[1][0], 'y': x[1][1], 'z': x[1][2], 'type':x[1][3] } ,  filtered_list )

		results, plot_stack = [], []

		# preparation for our calculations
		sin_theta = np.sin( theta )
		cos_theta = np.cos( theta )

		polygon_template = np.array(polygon_vertices)

		for step in range( steps ):
			distance_travelled = step * step_size

			position_x = x + distance_travelled * cos_theta
			position_y = y + distance_travelled * sin_theta 
			
			current_polygon_points =  polygon_template + [ position_x , position_y ]

			if return_plot_stack:
				# save the plot_stack:
				for i in range( len( current_polygon_points ) - 1 ):
				    vertex_set = [ current_polygon_points[i] ] + [ current_polygon_points[i+1] ]
				    plt_x, plt_y = zip(*vertex_set)
				    plot_stack.append([plt_x, plt_y, 'g'])

				vertex_set = [ current_polygon_points[0] ] + [ current_polygon_points[-1] ]
				plt_x, plt_y = zip(*vertex_set)
				plot_stack.append([plt_x, plt_y, 'g'])

			p = Path( current_polygon_points )

			# sample only those which are
			sample = filter( lambda r: p.contains_point( [ r['x'], r['y'] ] ) , r_vectors )
			results.append( ( distance_travelled , [ cell['id'] for cell in sample ] ) )


		# final points for plotting
		if return_plot_stack:
			plot_stack.append( [ [ x, x + step_size * steps * cos_theta], [ y, y + step_size * steps * sin_theta] ]  )

		if return_plot_stack:
			return results, plot_stack
		else:
			return results


	def cluster_search( self , *args, **kwargs):
		"""
			DEPRECATED PostProcess.sample_analyze()
		"""
		raise FutureWarning('This method is deprecated, use PostProcess.sample_analyze() instead')

	
	def sample_analyze( self , sample ):
		"""
			returns the 1D frequencies of each cluster 
			@params:
				sample: list of ( distance , [cell ids list] )
					this is often returned from a sampling strategy.
			@returns
				list of ( distance , frequency_analyze([ cell ids list] ) )
		"""
		return map( lambda cluster: ( cluster[0] , self.frequency_analyze( cluster[1] ) ) , sample )

	def frequency_analyze_ND ( self , clusters, **kwargs):
		"""
			gets the 2D frequency distribution between clusters (for now it works between two clusters only)
			@params:
				clusters / list of lists / [mandatory]
					a list of lists which contain the cellids to be analyzed. Examples:
						[ [ 1, 2, 3, 4, 5, 6 ] , [ 7, 8, 9, 10, 11, 12 ] ]
			@returns:
				2D np.array that is indexed by x,y = the frequency of mutations in respective clusters
				and entries representing the frequency of that count
		"""

		N = len( clusters )
		assert N == 2, 'frequency_analyze_ND does not support len(clusters) > 2 yet'

		freqs = []
		for cluster in clusters:
			freqs.append( self.frequency_analyze( cluster , return_loci = True )  )


		joint = [ ( freqs[0].get( k , 0 ) , freqs[1].get( k , 0 ) ) for k in set( freqs[0].keys() + freqs[1].keys() ) ] 
		self.joint = joint
		joint_frequency = Counter(joint)
		self.joint_frequency = joint_frequency # @TODO: remove this
		## first entry contains the most common element
		## second entry of the first entry contains the count of that element
		## those are the limits of our frequency diagram
		I = np.zeros( [ len( clusters[0] ) + 1 , len( clusters[1] ) + 1 ] , dtype=np.int_ )

		for k, v in joint_frequency.items():
			# reduce indicies by one to correct for the edges

			I[k[0]][k[1]] = v # frequency at that frequency

		return I, len( clusters[0] ) , len( clusters[1] )

	@staticmethod
	def plot_frequency_graph( frequency_results , title = '' , plot_stack = None , save_fig = None , scaled = False, ymax = None ):
		"""
			Plots the results of PostProcess.frequency_analyze()
			@params:
				frequency_results / results from PostProcess.frequency_analyze()
				title / str / ''
					title to append to the plot
				plot_stack / list 
					must contain arguments  to be executed for plotting by plt.plot
		"""

		frequency_results , number_of_cells = frequency_results

		# check if there are any results to show, otherwise just skip this.

		if len( frequency_results ):
			plt.figure()
			x, y = zip( *frequency_results.items() )
			x = np.array(x)
			x = x/float(number_of_cells)
			plt.plot(x, y, 'o')
			plt.ylabel('# of mutations shared')
			plt.xlabel('proportion of cells in cluster (total='+str(number_of_cells)+')')

			if title != '':
				title = '\n' + title

			plt.title('( # of shared mutations ) VS ( # of cells ) '+title)

			if scaled:
				if ymax is None:
					ymax = y.mean() + y.std()
				#end if

				plt.ylim( [ 0 , ymax ] )
				plt.xlim( [ 0 , 1 ] )

			#end if

			if plot_stack:
				while len( plot_stack ):
					args = plot_stack.pop()
					plt.plot(*args)
				#endwhile
			#endif

			if save_fig:
				plt.savefig( save_fig , format='png')
				plt.clf()
				plt.cla()
			else:
				plt.show()

	@staticmethod
	def plot_2D_frequency( frequency_results , title = '' , xlim = None , ylim = None , save_fig = None ):
		"""
			Plots the results of PostProcess.frequency_analyze_ND()
			@params:
				frequency_results / results from PostProcess.frequency_analyze_ND()
				title / str / ''
					title to append to the plot
		"""		

		frequency_results , N_1, N_2 = frequency_results

		# check if there are any results to show, otherwise just skip this
		if len( frequency_results ):
			# plt.imshow(  , interpolation = 'nearest' )
			plt.figure()
			
			plt.imshow( frequency_results , interpolation = 'nearest' , norm=LogNorm() , extent=[0, 1, 1, 0] )
			
			if xlim:
				plt.xlim( xlim )
			
			if ylim:
				plt.ylim( ylim )
			else:
				plt.ylim( tuple( reversed( plt.ylim() ) ) )


			cbar = plt.colorbar( orientation = 'vertical' )
			cbar.ax.set_ylabel('# of mutations')

			if title != '':
				title = '\n' + title

			plt.xlabel('proportion of cells in cluster 1 (N=' + str( N_1 ) + ')')
			plt.ylabel('proportion of cells in cluster 2 (N=' + str( N_2 ) +')')
			plt.title('2D Frequency Distribution of Inter-Cluster Mutations '+title)
			
			if save_fig:
				plt.savefig( save_fig , format='png')
				plt.clf()
				plt.cla()
			else:
				plt.show()


class TimeSeriesPlotters( object ):
	"""
		contains plotting functions for time series data
	"""

	@staticmethod
	def process_counts( cell_counts, dtype=np.int ):
		a = np.array( cell_counts, dtype = dtype )

		allcells = a[np.where(a[:,1] == -1 )]
		normal = a[np.where(a[:,1] == 0 )]
		cancer1 = a[np.where(a[:,1] == 1 )]
		cancer2 = a[np.where(a[:,1] == 2 )]

		t = normal[:,0]

		y0 = normal[:,2] 
		y1 = cancer1[:,2]
		y2 = cancer2[:,2]
		a = allcells[:,2]

		return t, y0, y1, y2, a


	@staticmethod
	def cell_counts_plot( cell_counts , file_name = './cell_count.png' ):		
		"""
			CELL NUMBER PLOTS
		"""

		plt.figure()

		t, y0, y1, y2, _ = TimeSeriesPlotters.process_counts( cell_counts )

		plt.plot(t, y0, label='normal cells')
		if not all(y1):
			plt.plot(t, y1, label='cancer1 cells')

		if not all(y2):
			plt.plot(t, y2, label='cancer2 cells')

		plt.xlabel('Time (MCS)')
		plt.ylabel('Number of cells')
		plt.title('Number of cells in tumor vs time')
		plt.legend()
		plt.savefig( file_name , format='png')
		plt.clf()
		plt.cla()

	@staticmethod
	def cell_proportions_plot( cell_counts , file_name = './cell_proportions.png' ):		
		"""
			CELL NUMBER PLOTS
		"""

		plt.figure()

		t, y0, y1, y2, a = TimeSeriesPlotters.process_counts( cell_counts )

		y0 = y0 / a.astype(float)
		y1 = y1 / a.astype(float)
		y2 = y2 / a.astype(float)

		plt.plot(t, y0, label='normal cells')
		if not all(y1):
			plt.plot(t, y1, label='cancer1 cells')

		if not all(y2):
			plt.plot(t, y2, label='cancer2 cells')

		plt.xlabel('Time (MCS)')
		plt.ylabel('Proportion of cells')
		plt.title('Proportion of cells in tumor vs time')

		plt.legend()
		plt.savefig( file_name , format='png')
		plt.clf()
		plt.cla()

	@staticmethod
	def all_cells_plot ( cell_counts , file_name = './all_cells_plot.png' ):

		t, y0, y1, y2, a = TimeSeriesPlotters.process_counts( cell_counts )

		plt.plot(t, a, )

		plt.xlabel('Time (MCS)')
		plt.ylabel('Number of cells')
		plt.title('Number of cells in tumor vs time')

		plt.savefig( file_name , format='png')
		plt.clf()
		plt.cla()



	@staticmethod
	def volume_plot( volume_data , file_name = './volume_plot.png' ):

		t, y0, y1, y2, a = TimeSeriesPlotters.process_counts( volume_data , dtype=np.float)


		plt.plot(t, y0, label='normal cells')
		if not all(y1):
			plt.plot(t, y1, label='cancer1 cells')

		if not all(y2):
			plt.plot(t, y2, label='cancer2 cells')

		plt.xlabel('Time (MCS)')
		plt.ylabel('Mean Volume of Cells')
		plt.title('Volume of cells in tumor vs time')
		plt.legend()
		plt.savefig( file_name , format='png')

		plt.clf()
		plt.cla()


def proportion_pairwise_differences ( allele_frequencies ):
	"""
		Calculates E(pi) = the expected proportion of pairwise
						differences 
		@params:
			allele_frequencies:
				the output from PostProcess.sample_analyze()
		@return:
			E_of_pi : calculated according to the equation:
						## E(pi) = sum_{l=1}^M 2 * X_l ( 1 - X_l)
						## l = locus , M = number of loci
						## X_l = pairwise proportional difference at site l

	"""
	frequency_results , number_of_cells = allele_frequencies

	if len( frequency_results ):
		# 
		x, _ = zip( *frequency_results.items() )
		x = np.array(x) 
		x = x / float(number_of_cells) # proportion of cells in cluster mutated.
		
		##
		## E(pi) = sum_{l=1}^M 2 * X_l ( 1 - X_l)
		## l = locus , M = number of loci
		## X_l = pairwise proportional difference at site l
		
		E_of_pi = np.sum( 2 * x * ( 1 - x ) ) 

		return E_of_pi 

	


def distance_between ( a , b ):
	q = b[0] - a[0] 
	p = b[1] - a[1] 
	r = b[2] - a[2]
	return math.sqrt( (p * p) + (q * q) + (r * r) )




