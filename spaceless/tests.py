from Model import Simulator

sim = Simulator()
# sim.run(time_steps=50)

timeseries = sim.run(time_steps=500,track_number=True)

import matplotlib.pyplot as plt

t, N = zip(*timeseries);

plt.figure()
plt.plot(t,N)
plt.show()


