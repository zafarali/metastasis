# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

import numpy as np
import matplotlib.pyplot as plt

# <codecell>

# generate some genomes
freq1 = np.random.randint(0, 500, size=10000).tolist()
freq2 = np.random.randint(0, 500, size=10000).tolist()
# genome1, genome2

# <codecell>

from collections import Counter

# <codecell>

# generate their counters:
g1_c = Counter(freq1)
g2_c = Counter(freq2)

# we already have these counts (i.e. counts of the # of mutations)

# <codecell>

results2 = [ ( g1_c.get( k , 0 ) , g2_c.get( k , 0 ) ) for k in set( g1_c.keys() + g2_c.keys() ) ] 
results2 = Counter(results2)

# <codecell>

#d = np.zeros([len(g1_c),len(g2_c)])
d = np.zeros([g1_c.most_common(1)[0][1],g2_c.most_common(1)[0][1]])
d.shape

# <codecell>

for k,v in results2.items():
    d[k[0]-1][k[1]-1] = v

# <codecell>

plt.imshow(d, interpolation='nearest')
plt.colorbar()
plt.show()

# <codecell>


# <codecell>


