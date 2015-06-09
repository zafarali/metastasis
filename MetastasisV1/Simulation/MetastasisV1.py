
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


from cc3dtools.Tracker import Tracker

tracker = Tracker(fileName='../division_output8615.csv')


from MetastasisV1Steppables import ConstraintInitializerSteppable
ConstraintInitializerSteppableInstance=ConstraintInitializerSteppable(sim,_frequency=1)
steppableRegistry.registerSteppable(ConstraintInitializerSteppableInstance)
        

from MetastasisV1Steppables import GrowthSteppable
GrowthSteppableInstance=GrowthSteppable(sim,_frequency=1)
steppableRegistry.registerSteppable(GrowthSteppableInstance)
        

from MetastasisV1Steppables import MitosisSteppable
MitosisSteppableInstance=MitosisSteppable(sim, tracker_instance=tracker, _frequency=10)
steppableRegistry.registerSteppable(MitosisSteppableInstance)
        

from MetastasisV1Steppables import DeathSteppable
DeathSteppableInstance=DeathSteppable(sim,_frequency=1)
steppableRegistry.registerSteppable(DeathSteppableInstance)

from MetastasisV1Steppables import ExtraMultiPlotSteppable
extraMultiPlotSteppable=ExtraMultiPlotSteppable(_simulator=sim,_frequency=10)
steppableRegistry.registerSteppable(extraMultiPlotSteppable)

CompuCellSetup.mainLoop(sim,simthread,steppableRegistry)
        
        