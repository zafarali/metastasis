import Post
from Model import Simulator, pDivisionFunction
import numpy as np

def build_SD(init_steps=200, post_steps=500, mean_mutations=50, \
	update_mean_mutations=120):

	
	sim = Simulator(mean_mutations=mean_mutations)

	sim.run(time_steps=init_steps, proportion_divide='auto_reduce', auto_reduce_magnitude=0.6)

	sim.create_cancer(update_mean_mutations=update_mean_mutations, p_division_function=pDivisionFunction.constant(0.5))

	sim.run(time_steps=500, stop_normal_divisions=True, proportion_divide='auto_reduce', auto_reduce_magnitude=0.6)

	sim.sort_genomes()

	return sim


def build_MD(init_steps=200, post_steps=500, mean_mutations=50, \
	update_mean_mutations=120):


	phenotype_template = {
	    'advantageous': (0, 0.5 * 100 )
	}

	sim = Simulator(mean_mutations=mean_mutations, phenotypes=phenotype_template)

	sim.run(time_steps=init_steps, proportion_divide='auto_reduce', auto_reduce_magnitude=0.6)

	sim.create_cancer(update_mean_mutations=update_mean_mutations)

	sim.run(time_steps=500, stop_normal_divisions=True, proportion_divide='auto_reduce', auto_reduce_magnitude=0.6)

	sim.sort_genomes()

	return sim

def build_CSC_reg(init_steps=200, post_steps=500, mean_mutations=50, \
	update_mean_mutations=120):
	

	# phenotype_template = {
	#     'advantageous': (0, 0.5 * 100 )
	# }

	sim = Simulator(mean_mutations=mean_mutations)

	sim.run(time_steps=init_steps, proportion_divide='auto_reduce', auto_reduce_magnitude=0.6)

	sim.create_CSC(update_mean_mutations=update_mean_mutations, p_division_function=pDivisionFunction.constant(0.5))

	sim.run(time_steps=500, stop_normal_divisions=True, proportion_divide='auto_reduce', auto_reduce_magnitude=0.6, age_mode=True)

	sim.sort_genomes()

	return sim



def regular_processor(sorted_sim, max_iteration = 10, iteration_magnitude=10, thresholds = [ 0.1, 0.5, 0.9, 1 ], subsample=0):
	stats = [ ('N', 't', 'S', 'SH', 'Epi', 'D') ]

	# sort genomes

	# obtain genomes
	g = sorted_sim.sorted_genomes

	for i in range(1,max_iteration):
		for t in thresholds:
			# repeated estimation
			to_be_averaged = []
			for r in range(0,20):
				to_be_processed = Post.split_genomes(g, N=iteration_magnitude*i, t=t)
				fa = Post.frequency_analyze(to_be_processed, subsample=subsample)
				to_be_averaged.append(Post.get_stats(*fa))

			avgd = tuple(np.mean(np.array(to_be_averaged), axis=0))

			this_iteration = (iteration_magnitude*i, t) + avgd
			stats.append(this_iteration)

	return stats





