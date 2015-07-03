## Author: @zafarali
## July 2015
# This is a pipeline for the analysis of the output files from CC3D simulations using the cc3dtools
# library.
# call as: python pipeline.py <analysis_directory>
# Analysis directory must have the following files:
# start_cells_*
# finish_cells_*
# division_*
# genomes_*
# specs.txt

import sys
import glob
import time
import os
import csv
import shutil
import json

event_sequence = 0

# convinience method to print things to the console with line numbers
def print2( string ):
	global event_sequence
	string = str(string)
	print ' Out [' + str(event_sequence) + ']: ' + string
	event_sequence += 1


time_info = '_'.join(time.asctime().split(' '))

# an analysis directory must be specified.
try:
	analysis_directory = './' + sys.argv[1]
	if not os.path.exists( analysis_directory ):
		raise IOError()
except IndexError:
	sys.exit('(!) No analysis_directory was specified: \nanalysis_directory must be the first argument')
except IOError:
	sys.exit('(!) Analysis directory does not exist!')

# read all necessary files
try:
	FILES = {
		'start' : glob.glob( analysis_directory + 'start_cells_*' )[0],
		'finish' : glob.glob( analysis_directory + 'finish_cells_*' )[0],
		'division' : glob.glob( analysis_directory + 'division_*' )[0],
		'genomes' : glob.glob( analysis_directory + 'genomes_*' )[0],
		'directory': analysis_directory,
		'specifications': glob.glob( analysis_directory + 'specs*' )[0],
		'out': analysis_directory+'/pipe_out_'+time_info
	}
	del analysis_directory
except IndexError:
	sys.exit('---> (!) Incomplete files for data analysis. Analysis directory must have the following files:\n\
	start_cells_*\n\
	finish_cells_*\n\
	division_*\n\
	genomes_*\n\
	specs.txt')


# load sampling strategies
print2('Importing Speficiations')

# will hold all specifications
SPECS = {}

with open( FILES['specifications'] , 'r' ) as f:
	reader = csv.reader( f )
	for row in reader:
		SPECS[ row[0] ] = SPECS.get( row[0] , {} )
		SPECS[ row[0] ][ row[1] ] = row[2:]

# convinience method for looking up SPEC object
def SPEC_lookup( function_name , argument_name , default = None ):
	fn_specs = SPECS.get( function_name , None )
	if fn_specs is None:
		return default
	else:
		return fn_specs.get( argument_name , default )


print2('SPECS LOADED')

if SPEC_lookup( 'global', 'sample' ):
	if os.path.exists( FILES['directory'] + 'sampling.json' ):
		FILES['sampling'] = FILES['directory'] + 'sampling.json'

		with open( FILES['directory'] + 'sampling.json' ) as f:
			sampling_strategies = json.load( f )['sampling_strategies']
			print2('sampling strategies loaded')
	else:
		sys.exit('--->(!) No sampling strategy was provided in the form of sampling.json')


# create the output directory and save specs and sampling there
if not os.path.exists( FILES['out'] ):
    os.makedirs( FILES['out'] )

shutil.copy2( FILES['specifications'] , FILES['out'] + '/specs.txt' )

if SPEC_lookup('global', 'sample'):
	shutil.copy2( FILES['sampling'] , FILES['out'] + '/sampling.json' )

print2('File paths saved, directory created.')

print2('saving to:' + FILES['out'])



from cc3dtools.PostProcess import SpacePlot, PostProcess
from cc3dtools.GenomeCompare import GenomeCompare
from cc3dtools.Lineage import MultiLineage, lineage_bfs

sp = SpacePlot( start_file = FILES['start'], \
	end_file = FILES['finish'], \
	projection = SPEC_lookup( 'spatial_plot' , 'projection', '2d' ) \
	)


"""
S P A C E P L O T
"""
if SPEC_lookup('global', 'spatial_plot'):
	sp.plot_all( hide_numbers = SPEC_lookup( 'spatial_plot', 'hide_numbers', False ) , \
		save_fig = FILES['out']+'/spatial_plot.png' \
		)
	print2('saved spatial_plot')


"""
L I N E A G E S
"""
ml = MultiLineage( file_name = FILES['division'] )

if SPEC_lookup( 'global' , 'multilineage' ):
	if SPEC_lookup('multilineage', 'color_by_branch') or SPEC_lookup('multilineage', 'save_all'):
		for i , lineage in enumerate(ml.lineages):
			current_dir = FILES['out'] + '/detail_lineage' + str( i + 1 )
			os.makedirs( current_dir )
			levels = lineage_bfs( lineage['lineage'] )

			if SPEC_lookup('multilineage', 'color_by_branch'):
				for depth, level in levels.items():
					sp.plot_selected( *level , depth = depth , save_fig= current_dir + '/depth'+str(depth)+'.png' )
					print2('saved color_by_branch for lineage: ' + str(i+1) + 'depth:' + str(depth) )

			if SPEC_lookup('multilineage', 'save_all'):
				lineage['lineage'].draw( width=20, save_fig = current_dir + '/basic.png' )
				print2('saved overall lineage for '+str(i+1))

	if SPEC_lookup('multilineage', 'color_by_initial'):
		member_lists = map( lambda lineage: lineage['members'], ml.lineages )
		# print member_lists[0]
		sp.plot_selected( *member_lists, depth=1, save_fig = FILES['out'] + '/space_plot_by_initial.png' )
		print2('saved space_plot_by_initial')

# if SPEC_lookup( 'global' , 'sample' )




