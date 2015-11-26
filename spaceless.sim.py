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
		'mean_mutations': 1,
		'update_mean_mutations':25,
		'init_steps': INIT_STEPS,
		'post_steps': POST_STEPS,
		'auto_reduce_mangitude': AUTO_REDUCE_MAGNITUDE
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

	file_name = './spaceless_data/'+MODEL_TYPE+'.'+INIT_STEPS+'.'+POST_STEPS+'.'+AUTO_REDUCE_MAGNITUDE+'.'+time.ctime()

	vals = regular_processor(sim, subsample=100, thresholds=[0.7,0.8,0.9,1])
	print sim.sorted_genomes
	print 'total cells:',len(sim.cells)
	print vals

	with open(file_name+'stats.csv', 'w') as f:
		writer = csv.writer(f)
		for val in vals:
			writer.writerow(val)



