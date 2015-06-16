## Author @zafarali
## June 2015

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import csv
import math
import time
import pickle


def spatial_plot( start_file = None , end_file = None , type_colors = ( 'r', 'b', 'g' ) , format = ( 'id' , 'type' , 'x', 'y', 'z' ) , projection = '2d' ):
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
				ax.text( data['x']+0.5 , data['y']+0.5 , data['z'] , str( data['id'] ) , horizontalalignment = 'center' , color = type_colors[ data['type'] - 1 ] )
			else:
				plt.plot( data['x'] , data['y'] , marker+type_colors[ data['type'] - 1 ] )
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
			(!) warning, this is a computationally heavy function

		"""
		assert self.__executed__ == False, 'Cannot re-execute an already executed PostProcess object'
		pass

		results = {}
		
		num_loops = 0

		start_time = time.time()
		print 'start time:',start_time
		for index, i in enumerate(self.cells_available):

			results[i] = {}

			for j in self.cells_available:
				
				num_loops += 1

				if i == j: # do not compute when i == j as they are just 0
					shared = len( self.gc.genomes[i].get_mutated_loci() )
					private = 0
					distance = 0
					results[i][j] = { 'distance' : distance , 'private' : private , 'shared' : shared }

					# we reach the middle of the matrix, thus we break and skip to the next coloumn.
					# print 'completed ',i
					break

				a = self.cell_locations[i]
				b = self.cell_locations[j]

				distance = distance_between( a , b )
				shared = self.gc.comm( i , j )
				private = self.gc.diff( i , j ) 

				# the resulting matrix is going to be symmetric
				results[i][j] = { 'distance' : distance , 'private' : private , 'shared' : shared }
				results[j][i] = { 'distance' : distance , 'private' : private , 'shared' : shared }

			# end for j
		# end for i

		self.data = results
		self.__executed__ = True

		end_time = time.time()
		
		print 'Completed\n Total Number of loops:' + str( num_loops )
		print 'Total Time:' + str( end_time - start_time ) + 's'

	def pwise ( self , c1 , c2 ):
		"""
			Returns pairwise information for cells c1 and c2
		"""
		assert self.__executed__ == True, 'You must first PostProcess.execute() before you can access other methods'

		return self.data[c1][c2]

		pass
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




