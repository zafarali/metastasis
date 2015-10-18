
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




from LimboMutationResponseSteppables import ConstraintInitializerSteppable
ConstraintInitializerSteppableInstance=ConstraintInitializerSteppable(sim,_frequency=1)
steppableRegistry.registerSteppable(ConstraintInitializerSteppableInstance)
        

from LimboMutationResponseSteppables import GrowthSteppable
GrowthSteppableInstance=GrowthSteppable(sim,_frequency=1)
steppableRegistry.registerSteppable(GrowthSteppableInstance)


from LimboMutationResponseSteppables import UtillitySteppable
UtillitySteppableInstance=UtillitySteppable(sim,_frequency=1)
steppableRegistry.registerSteppable(UtillitySteppableInstance)


from LimboMutationResponseSteppables import MitosisSteppable
MitosisSteppableInstance=MitosisSteppable(sim, _frequency=10)
steppableRegistry.registerSteppable(MitosisSteppableInstance)


from LimboMutationResponseSteppables import SuperTracker
SuperTrackerInstance = SuperTracker(sim, _frequency= 100 )
steppableRegistry.registerSteppable(SuperTrackerInstance)


from LimboMutationResponseSteppables import DeathCheckSteppable
DeathCheckSteppableInstance = DeathCheckSteppable(sim, _frequency= 10 )
steppableRegistry.registerSteppable(DeathCheckSteppableInstance)


from LimboMutationResponseSteppables import DeathSteppable
DeathSteppableInstance = DeathSteppable(sim, _frequency= 1 )
steppableRegistry.registerSteppable(DeathSteppableInstance)



CompuCellSetup.mainLoop(sim,simthread,steppableRegistry)
        
        