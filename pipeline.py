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


time_info = '_'.join(time.asctime().split(' '))
try:
	analysis_directory = './' + sys.argv[1]
except IndexError:
	sys.exit('--->No analysis_directory was specified: \nanalysis_directory must be the first argument')

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

if not os.path.exists( FILES['out'] ):
    os.makedirs( FILES['out'] )

print 'Importing specifications'


# will hold all specifications
SPECS = {}

with open( FILES['specifications'] , 'r' ) as f:
	reader = csv.reader( f )
	for row in reader:
		SPECS[ row[0] ] = SPECS.get( row[0] , {} )
		SPECS[ row[0] ][ row[1] ] = bool( row[2] )

# convinience method for looking up SPEC object
def SPEC_lookup( function_name , argument_name , default = None ):
	fn_specs = SPECS.get( function_name , None )
	if fn_specs is None:
		return default
	else:
		return fn_specs.get( argument_name , default )

shutil.copy2( FILES['specifications'] , FILES['out'] + '/specs.txt' )

print 'SPECS LOADED:'
print SPECS

print 'File paths saved, directory created.'
print 'saving to:' + FILES['out']



# spatial_plot( start_file = FILES['start'], \
# 	end_file = FILES['finish'], \
# 	hide_numbers = True, \
# 	save_fig = FILES['directory']+'/test_save.png'
# 	)