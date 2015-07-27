#!/usr/bin/env python
import sys
sys.path.append('/Users/zafaraliahmed/summer15/metastasis/')
import dispatcher

cpu=int(sys.argv[2]) #number of cpus to use


# stores the commands we are going to run
commandlist=[] 


# names of the simulations we are going to run
simulation_names = sys.argv[3:]

# number of times each simulation must be run
num_times = int( sys.argv[1] )

for simulation in simulation_names: 
	for i in range(num_times):
		#for bootnum in range(100):
		commandlist.append('~/CC3D_3.7.3/runScript.command -i ~/summer15/metastasis/'+simulation+'/'+simulation+'.cc3d')
		
dispatcher.dispatcher(commandlist,slots=min(len(commandlist),cpu))