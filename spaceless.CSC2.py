print """Cancer Stem Cell with driver support Simulation
Spaceless Model
"""

from spaceless.Toys import build_CSC_reg, regular_processor
from spaceless import Post
from spaceless.Model import pDivisionFunction
import time
import matplotlib.pyplot as plt
import numpy as np
import csv

sim = build_CSC_reg(mean_mutations=1, init_steps=500, update_mean_mutations=25, p_division_function=pDivisionFunction.sigmoid(), post_steps=500)

# from cc3dtools.Genome import save_genomes2, get2_to_dict
# from cc3dtools.GenomeCompare import GenomeCompare


file_name = './spaceless_data/CSC_MD.'+time.ctime()

# # # genomes =  sim.get_genomes()
# # types = sim.get_types()
# # sim.sort_genomes()
# # genomes = sim.sorted_genomes
print sim.cell_stats()

vals = regular_processor(sim, subsample=100)

# print vals


# vals = regular_processor(sim, thresholds=[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1])

statistics = np.array(vals[1:])
plt.figure('DvsN')
to_plot = 5
threshold = np.where(statistics[:,1] == 0.1)
plt.scatter(statistics[threshold,0], statistics[threshold,to_plot], label='0.1', c='r')
threshold = np.where(statistics[:,1] == 0.5)
plt.scatter(statistics[threshold,0], statistics[threshold,to_plot], label='0.5', c='g')
threshold = np.where(statistics[:,1] == 0.9)
plt.scatter(statistics[threshold,0], statistics[threshold,to_plot], label='0.9', c='b')
threshold = np.where(statistics[:,1] == 1)
plt.scatter(statistics[threshold,0], statistics[threshold,to_plot], label='1', c='y')
plt.title('D vs sample size (for different thresholds')
plt.legend()
plt.show()

plt.figure('EpivsN')
to_plot = 4
threshold = np.where(statistics[:,1] == 0.1)
plt.scatter(statistics[threshold,0], statistics[threshold,to_plot], label='0.1', c='r')
threshold = np.where(statistics[:,1] == 0.5)
plt.scatter(statistics[threshold,0], statistics[threshold,to_plot], label='0.5', c='g')
threshold = np.where(statistics[:,1] == 0.9)
plt.scatter(statistics[threshold,0], statistics[threshold,to_plot], label='0.9', c='b')
threshold = np.where(statistics[:,1] == 1)
plt.scatter(statistics[threshold,0], statistics[threshold,to_plot], label='1', c='y')
plt.title('Epi vs sample size (for different thresholds')
plt.legend()
plt.show()

plt.figure('SHvsN')
to_plot = 3
threshold = np.where(statistics[:,1] == 0.1)
plt.scatter(statistics[threshold,0], statistics[threshold,to_plot], label='0.1', c='r')
threshold = np.where(statistics[:,1] == 0.5)
plt.scatter(statistics[threshold,0], statistics[threshold,to_plot], label='0.5', c='g')
threshold = np.where(statistics[:,1] == 0.9)
plt.scatter(statistics[threshold,0], statistics[threshold,to_plot], label='0.9', c='b')
threshold = np.where(statistics[:,1] == 1)
plt.scatter(statistics[threshold,0], statistics[threshold,to_plot], label='1', c='y')
plt.title('SH vs sample size (for different thresholds')
plt.legend()
plt.show()


plt.figure('Dvst')
to_plot = 5
threshold = np.where(statistics[:,0] == 525)
plt.scatter(statistics[threshold,1], statistics[threshold,to_plot], label='N=525', c='r')
threshold = np.where(statistics[:,0] == 775)
plt.scatter(statistics[threshold,1], statistics[threshold,to_plot], label='N=775', c='g')
threshold = np.where(statistics[:,0] == 1025)
plt.scatter(statistics[threshold,1], statistics[threshold,to_plot], label='N=1025', c='b')
threshold = np.where(statistics[:,0] == 1900)
plt.scatter(statistics[threshold,1], statistics[threshold,to_plot], label='N=1900', c='y')
plt.title('D vs t (for different N)')
plt.legend()
plt.show()

plt.figure('Svst')
to_plot = 3
threshold = np.where(statistics[:,0] == 525)
plt.scatter(statistics[threshold,1], statistics[threshold,to_plot], label='N=525', c='r')
threshold = np.where(statistics[:,0] == 775)
plt.scatter(statistics[threshold,1], statistics[threshold,to_plot], label='N=775', c='g')
threshold = np.where(statistics[:,0] == 1025)
plt.scatter(statistics[threshold,1], statistics[threshold,to_plot], label='N=1025', c='b')
threshold = np.where(statistics[:,0] == 1900)
plt.scatter(statistics[threshold,1], statistics[threshold,to_plot], label='N=1900', c='y')
plt.title('SH vs t (for different N)')
plt.legend()
plt.show()

plt.figure('epivst')
to_plot = 4
threshold = np.where(statistics[:,0] == 525)
plt.scatter(statistics[threshold,1], statistics[threshold,to_plot], label='N=525', c='r')
threshold = np.where(statistics[:,0] == 775)
plt.scatter(statistics[threshold,1], statistics[threshold,to_plot], label='N=775', c='g')
threshold = np.where(statistics[:,0] == 1025)
plt.scatter(statistics[threshold,1], statistics[threshold,to_plot], label='N=1025', c='b')
threshold = np.where(statistics[:,0] == 1900)
plt.scatter(statistics[threshold,1], statistics[threshold,to_plot], label='N=1900', c='y')
plt.title('Pi vs t (for different N)')
plt.legend()
plt.show()


# with open(file_name+'stats.csv', 'w') as f:
# 	writer = csv.writer(f)
# 	for val in vals:
# 		writer.writerow(val)


# for testing purposes
# this should show that any cell with max division of 4
# should not have the middle element larger than 4
# from collections import Counter
# print Counter(map(lambda cell: (cell.max_divisions, cell.number_of_divisions, cell.cell_type), sim.cells.values()))