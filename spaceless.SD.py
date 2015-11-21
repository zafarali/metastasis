print """Single Driver Simulation
Spaceless Model
"""

from spaceless.Toys import build_SD, regular_processor
from spaceless import Post
import time
import matplotlib.pyplot as plt
import numpy as np
import csv

sim = build_SD(mean_mutations=1, update_mean_mutations=25, post_steps=1500)

# from cc3dtools.Genome import save_genomes2, get2_to_dict
# from cc3dtools.GenomeCompare import GenomeCompare



file_name = './spaceless_data/SD.'+time.ctime()

# # genomes =  sim.get_genomes()
# types = sim.get_types()
# sim.sort_genomes()
# genomes = sim.sorted_genomes


vals = regular_processor(sim)

with open(file_name+'stats.csv', 'w') as f:
	writer = csv.writer(f)
	for val in vals:
		writer.writerow(val)

