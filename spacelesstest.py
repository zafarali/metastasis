import numpy as np

np.random.seed(123)

from spaceless.Model import Simulator


phenotype_template = {
    'advantageous': (0, 0.5 * 100 )
}

sim = Simulator( phenotypes=phenotype_template, cell_attributes = { 'mean_mutations' : 1, 'chromosome_order':2 , 'ploidy':2} )

sim.run(time_steps=10)

cancer_id = sim.create_cancer(update_mean_mutations=5)

print 'cancer created: ',sim.cells[cancer_id].is_cancer()
print 'cancer id: ',cancer_id

# sim.run(time_steps=100, stop_normal_divisions=True, proportion_divide='auto_reduce', auto_reduce_magnitude=0.75)

# print sim.cells
# print sim.cells.items()[-1][1].cell_type

from spaceless.Model import SelectionDistribution


print np.unique( SelectionDistribution.equal( sim.cells.values() ) )
# print map( lambda cell: cell.p_division_function, sim.cells.values() )

print np.unique( SelectionDistribution.cancer_only( sim.cells.values() ) )
selected_cell = sim.cells.values()[-1]
print selected_cell.phenotype.get_counts()['advantageous']
print selected_cell.p_division()

print sim.cells.values()[0].p_division()
print sim.cells.values()[cancer_id].p_division()
# print map()

# print selected_cell.genome.get_mutated_loci()
print selected_cell.genome.mutation_rate

from cc3dtools.Genome import save_genomes2

save_genomes2(map(lambda cell: cell.genome, sim.cells.values()))