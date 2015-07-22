
import sys
from os import environ
from os import getcwd
import string
import getpass

sys.path.append(environ["PYTHON_MODULE_PATH"])

user = getpass.getuser()
sys.path.append('/Users/'+user+'/summer15/metastasis/')


import CompuCellSetup


sim,simthread = CompuCellSetup.getCoreSimulationObjects()
            
# add extra attributes here
            
CompuCellSetup.initializeSimulationObjects(sim,simthread)
# Definitions of additional Python-managed fields go here
        
#Add Python steppables here
steppableRegistry=CompuCellSetup.getSteppableRegistry()

# holder = {}




from LimboSteppables import ConstraintInitializerSteppable
ConstraintInitializerSteppableInstance=ConstraintInitializerSteppable(sim,_frequency=1)
steppableRegistry.registerSteppable(ConstraintInitializerSteppableInstance)
        

from LimboSteppables import GrowthSteppable
GrowthSteppableInstance=GrowthSteppable(sim,_frequency=1)
steppableRegistry.registerSteppable(GrowthSteppableInstance)
        

from LimboSteppables import MitosisSteppable
MitosisSteppableInstance=MitosisSteppable(sim, _frequency=10)
steppableRegistry.registerSteppable(MitosisSteppableInstance)
        


# from LimboSteppables import DeathSteppable, DeathCheckSteppable

# DeathSteppableInstance=DeathSteppable(sim, _frequency=1)
# DeathCheckSteppableInstance=DeathCheckSteppable(sim, _frequency=10)

# steppableRegistry.registerSteppable(DeathSteppableInstance)
# steppableRegistry.registerSteppable(DeathCheckSteppableInstance)

from LimboSteppables import SuperTracker
SuperTrackerInstance = SuperTracker(sim, _frequency= 10 )
steppableRegistry.registerSteppable(SuperTrackerInstance)


# from LimboSteppables import ExtraMultiPlotSteppable
# extraMultiPlotSteppable=ExtraMultiPlotSteppable(_simulator=sim,_frequency=10)
# steppableRegistry.registerSteppable(extraMultiPlotSteppable)


CompuCellSetup.mainLoop(sim,simthread,steppableRegistry)
        
        