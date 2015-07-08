"""
Example Usage
#!/usr/bin/env python
import sys
sys.path.append("/Users/gravel/bin/scripts/dispatcher/")
import dispatcher
cpu=6 #number of cpus to use
commandlist=[]
 

for cohort in ["SCCS","HRS_current"]: 
	for model in ["pp_fix","pp_px_fix","ppx_xxp_fix","ppx_xxp_xxp_fix"]:
		#for bootnum in range(100):
		commandlist.append("./run_tracts_boot.py "+model+" "+cohort)
		
dispatcher.dispatcher(commandlist,slots=min(len(commandlist),cpu))

"""
import subprocess
from time import sleep
import numpy

def dispatcher(comm,slots=4,delay=2,nice=None,verbose=100):
	"""launches processes in list "comm" using the subprocess module. 
	Makes sure no more than "slots" processes are launched at any given time. 
	Waits a time "delay" seconds between checks of process completion. 
	Returns a list containing the stdout of each process. 
	Verbose: approximately "verbose" reports""" 
	
	
	nProc=len(comm)
	ploteach=(nProc/verbose)+1
	print "will outpout every %d iteration" % (ploteach,) 
	readout=[]
	a=[]
	for it in range(min(slots,nProc)):
		a.append(0)
	for it in range(0,min(slots,nProc)):
		
		if nice is not None:
			comm[it]='nice -%d' % nice[it]
		
		a[it]=(it,subprocess.Popen(comm[it],shell=True,stdout=subprocess.PIPE))
	
	for it in range(slots,nProc):	
		done=False
		while done is False:
			for slotNo in range(0,slots):
				if a[slotNo][1].poll() is not None:
					place=slotNo
					done=True
					break
					sleep(delay)
			sleep(1)#To make sure dispatcher is not taking a cpu by itself (may be the case with certain version of python.
			
		if done==True:
				
			readout.append((a[place][0],a[place][1].stdout.readlines()))
			if nice is not None:
				comm[it]='nice -%d' % nice[place]	
			if it % ploteach==0:
				print "sending job"
				print comm[it]
			a[place]=(it,subprocess.Popen(comm[it],shell=True,stdout=subprocess.PIPE))
				
	#wait for the last slots processes
	for process in a:
		print process
		process[1].wait()
		print "finished job ", process[1]
		readout.append((process[0],process[1].stdout.readlines()))
		
			
	return readout