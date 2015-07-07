
import sys
from os import environ
from os import getcwd
import string

sys.path.append(environ["PYTHON_MODULE_PATH"])
sys.path.append('/Users/zafaraliahmed/summer15/metastasis/')


import CompuCellSetup


sim,simthread = CompuCellSetup.getCoreSimulationObjects()
            
# add extra attributes here
            
CompuCellSetup.initializeSimulationObjects(sim,simthread)
# Definitions of additional Python-managed fields go here
        
#Add Python steppables here
steppableRegistry=CompuCellSetup.getSteppableRegistry()

# holder = {}




from MetastasisMutationResponseSteppables import ConstraintInitializerSteppable
ConstraintInitializerSteppableInstance=ConstraintInitializerSteppable(sim,_frequency=1)
steppableRegistry.registerSteppable(ConstraintInitializerSteppableInstance)
        

from MetastasisMutationResponseSteppables import GrowthSteppable
GrowthSteppableInstance=GrowthSteppable(sim,_frequency=1)
steppableRegistry.registerSteppable(GrowthSteppableInstance)
        

from MetastasisMutationResponseSteppables import MitosisSteppable
MitosisSteppableInstance=MitosisSteppable(sim, _frequency=10)
steppableRegistry.registerSteppable(MitosisSteppableInstance)
        

from MetastasisMutationResponseSteppables import DeathSteppable
DeathSteppableInstance=DeathSteppable(sim,_frequency=1)
steppableRegistry.registerSteppable(DeathSteppableInstance)

# 

from MetastasisMutationResponseSteppables import NeighborTrackerPrinterSteppable
steppableRegistry.registerSteppable( NeighborTrackerPrinterSteppable(_simulator = sim, _frequency = 10 ) )


"""
	COMMENT OUT DUE TO CLI MODE
"""
# from MetastasisV1Steppables import ExtraMultiPlotSteppable
# extraMultiPlotSteppable=ExtraMultiPlotSteppable(_simulator=sim,_frequency=10)
# steppableRegistry.registerSteppable(extraMultiPlotSteppable)

CompuCellSetup.mainLoop(sim,simthread,steppableRegistry)
        
        