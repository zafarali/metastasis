## Author @zafarali
## June 2015

import matplotlib.pyplot as plt
import csv

def spatial_plot( start_file = None, end_file = None, type_colors = ( 'r', 'b', 'g' ), format = ( 'id' , 'type' , 'x', 'y', 'z' ) ):
	assert start_file is not None and end_file is not None, 'must specify files'

	points = []

	starter_cells = [] # stores as (cellid, type)
	with open( start_file , 'r' ) as f:
		reader = csv.reader( f )
		for row in reader:
			starter_cells.append(  int( row[0] )  )

	plt.figure()

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
			print data

			plt.text( data['x'] , data['y'] , data['id'] , horizontalalignment = 'center', color = type_colors[ data['type'] - 1 ] )

	plt.show()







