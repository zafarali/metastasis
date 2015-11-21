import Post
from Model import Simulator, pDivisionFunction

def build_SD(init_steps=200, post_steps=500, mean_mutations=50, \
	update_mean_mutations=120):

	
	sim = Simulator(mean_mutations=mean_mutations)

	sim.run(time_steps=init_steps, proportion_divide='auto_reduce', auto_reduce_magnitude=0.6)

	sim.create_cancer(update_mean_mutations=update_mean_mutations, p_division_function=pDivisionFunction.constant(0.5))

	sim.run(time_steps=500, stop_normal_divisions=True, proportion_divide='auto_reduce', auto_reduce_magnitude=0.6)


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

	return sim


def regular_processor(sim, max_iteration = 10, iteration_magnitude=10, thresholds = [ 0.1,0.25,0.5,0.75,0.9,1 ]):
	stats = [ ('N', 't', 'S', 'SH', 'D', 'Epi') ]

	sim.sort_genomes()

	g = sim.sorted_genomes

	for i in range(0,25):
		for t in thresholds:
			# repeated estimation
			for r in range(0,20):
				Post.