#!/usr/bin/env python

from spaceless import Toys
from spaceless.Model import pDivisionFunction
from spaceless import Post
import time
import matplotlib.pyplot as plt
import numpy as np
import csv
import sys

if len(sys.argv) < 6:
	print """
		You must supply the following:
		MODEL_TYPE: {CSCSD, MD, SD, CSCMD}
		INIT_STEPS: INT
		POST_STEPS: INT
		AUTO_REDUCE_MAGNITUDE
		SAVE_DIR: DIRECTOR
	"""
else:

	MODEL_TYPE = str(sys.argv[1])
	INIT_STEPS = int(sys.argv[2])
	POST_STEPS = int(sys.argv[3])
	AUTO_REDUCE_MAGNITUDE = float(sys.argv[4])
	SAVE_DIR = str(sys.argv[5])

	arguments = {
		'mean_mutations': 50,
		'update_mean_mutations':120,
		'init_steps': INIT_STEPS,
		'post_steps': POST_STEPS,
		'auto_reduce_magnitude': AUTO_REDUCE_MAGNITUDE
	}
	print 'MODEL_TYPE',MODEL_TYPE
	print 'ARGUMENTS',arguments

	if MODEL_TYPE == 'SD':
		builder = Toys.build_SD
	elif MODEL_TYPE == 'MD':
		builder = Toys.build_MD
	elif MODEL_TYPE == 'CSCMD':
		builder = Toys.build_CSC_reg
		arguments.update({
			'p_division_function': pDivisionFunction.sigmoid()
			})
	elif MODEL_TYPE == 'CSCSD':
		builder = Toys.build_CSC_reg

	sim = builder(**arguments)
	sim.sort_genomes()

	file_name = './spaceless_data/'+MODEL_TYPE+'.'+str(INIT_STEPS)+'.'+str(POST_STEPS)+'.'+str(AUTO_REDUCE_MAGNITUDE)+'.'+time.ctime()

	vals = Toys.regular_processor(sim, subsample=100, max_iteration=150, iteration_magnitude=100, thresholds=[0, 0.00001])
	# print sim.sorted_genomes
	print 'total cells:',len(sim.cells)
	print '-----'
	print vals

	with open(file_name+'stats.csv', 'w') as f:
		writer = csv.writer(f)
		for val in vals:
			writer.writerow(val)



