#!/usr/bin/env python
import sys
sys.path.append('/Users/zafaraliahmed/summer15/metastasis/')
import dispatcher

cpu=4 #number of cpus to use


# stores the commands we are going to run
commandlist=[] 


# names of the simulations we are going to run
simulation_directories = sys.argv[1:]



for simulation_directory in simulation_directories: 
	commandlist.append('python pipeline.py ./simulation_out/'+simulation_directory+'/')
		
dispatcher.dispatcher(commandlist,slots=min(len(commandlist),cpu))