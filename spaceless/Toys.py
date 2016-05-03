import Post
from Model import Simulator, pDivisionFunction
import numpy as np
import time
import sys

def build_SD(init_steps=200, post_steps=500, mean_mutations=50, \
	update_mean_mutations=120, auto_reduce_magnitude=0.6, chromosome_order=15):
	# phenotype_template = {
	#     'advantageous': (0, 0.5 * 10**chromosome_order )
	# }
	sim = Simulator(mean_mutations=mean_mutations, chromosome_order=chromosome_order)

	sim.run(time_steps=init_steps, proportion_divide='auto_reduce', auto_reduce_magnitude=auto_reduce_magnitude)

	sim.create_cancer(update_mean_mutations=update_mean_mutations, p_division_function=pDivisionFunction.constant(0.5))

	sim.run(time_steps=post_steps, stop_normal_divisions=True, proportion_divide='auto_reduce', auto_reduce_magnitude=auto_reduce_magnitude)

	sim.sort_genomes()

	return sim


def build_MD(init_steps=200, post_steps=500, mean_mutations=50, \
	update_mean_mutations=120, auto_reduce_magnitude=0.6, chromosome_order=15):


	phenotype_template = {
	    'advantageous': (0, 0.1 * 10**chromosome_order )
	}

	sim = Simulator(mean_mutations=mean_mutations, phenotypes=phenotype_template, chromosome_order=chromosome_order)

	sim.run(time_steps=init_steps, proportion_divide='auto_reduce', auto_reduce_magnitude=auto_reduce_magnitude)

	sim.create_cancer(update_mean_mutations=update_mean_mutations)

	sim.run(time_steps=post_steps, stop_normal_divisions=True, proportion_divide='auto_reduce', auto_reduce_magnitude=auto_reduce_magnitude)

	sim.sort_genomes()

	return sim

def build_CSC_reg(init_steps=200, post_steps=500, mean_mutations=50, \
	update_mean_mutations=120, p_division_function=pDivisionFunction.constant(0.5), \
	auto_reduce_magnitude=0.6, chromosome_order=15):
	

	phenotype_template = {
	    'advantageous': (0, 0.1 * 10**chromosome_order )
	}

	sim = Simulator(mean_mutations=mean_mutations, phenotypes=phenotype_template, chromosome_order=chromosome_order)

	sim.run(time_steps=init_steps, proportion_divide='auto_reduce', auto_reduce_magnitude=auto_reduce_magnitude)

	sim.create_CSC(update_mean_mutations=update_mean_mutations, p_division_function=p_division_function)

	sim.run(time_steps=post_steps, stop_normal_divisions=True, proportion_divide='auto_reduce', auto_reduce_magnitude=auto_reduce_magnitude)

	sim.sort_genomes()

	return sim



def regular_processor(sorted_sim, max_iteration = 100, iteration_magnitude=25, thresholds = [ 0.1, 0.5, 0.9, 1 ], subsample=0):
	stats = [ ('N', 't', 'S', 'SH', 'Epi', 'D', 'proportion_cancer', 'sd_S', 'sd_SH', 'sd_Epi', 'sd_D', 'sd_in_proportion' ) ]

	# sort genomes

	# obtain genomes
	g = sorted_sim.sorted_genomes
	print('Processing started at '+str(time.ctime()))
	start_time = time.time()
	for i in range(1,max_iteration,5):
		print 'Sampling N=',iteration_magnitude*i
		for t in thresholds:
			# repeated estimation
			to_be_averaged = []
			for r in range(0,10):
				to_be_processed, N_real, proportion_cancer = Post.split_genomes(g, N=iteration_magnitude*i, t=t)
				fa = Post.frequency_analyze(to_be_processed, subsample=subsample)
				to_be_averaged.append( Post.get_stats(*fa) + ( proportion_cancer, ) )
				# to_be_averaged.append(proportion_cancer)
			print '->Completed processing a threshold=',t
			avgd = tuple(np.mean(np.array(to_be_averaged), axis=0))
			sds = tuple(np.std(np.array(to_be_averaged), axis=0))

			this_iteration = (N_real, t) + avgd + sds
			stats.append(this_iteration)
		print 'Completed sampling of N=',iteration_magnitude*i
		sys.stdout.flush()
	print('regular_processor successfully completed at '+str(time.ctime()))
	print('Total time taken:'+str((time.time() - start_time))+'seconds')
	return stats





