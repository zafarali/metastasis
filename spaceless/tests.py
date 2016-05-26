
from Model import Chromosome

c = Chromosome(mean_mutations=5)

c.mutate()
c.mutate()

mutated_loci = list(c.get_mutated_loci())


to_check = mutated_loci[0]

to_check_x = mutated_loci[0].to_int()


assert c.is_mutated(to_check) == True, 'Mutation() equality failed'
assert c.is_mutated(to_check_x) == True, 'integer loci equality failed'