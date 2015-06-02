
from PySteppables import *
import CompuCell
import sys

import numpy as np

from PySteppablesExamples import MitosisSteppableBase
            
# phenotypes = {
#     'deleterious'
# }

divide_times = {'last_division':0}

class ConstraintInitializerSteppable(SteppableBasePy):
    def __init__(self,_simulator,_frequency=1):
        SteppableBasePy.__init__(self,_simulator,_frequency)
    def start(self):
        print type(self.cellList)
        print dir(self.cellList)
        print type(self.cellList.inventory)
        print dir(self.cellList.inventory)
        num_cells = len(self.cellList)


        for cell in self.cellList:
            
            divide_times[cell.id] = 0

            cell.targetVolume=25
            cell.lambdaVolume=2.0
            # holder[cell.id] = { 'g': Genome( mutation_rate = 20 , genome_order = 10 ), 'p': Phenotype() }

        r = np.random.randint(0,num_cells)
        self.cellList.inventory.attemptFetchingCellById(r).type = self.CANCER1
        self.cellList.inventory.attemptFetchingCellById(r).targetVolume = 35
        

        

class GrowthSteppable(SteppableBasePy):
    def __init__(self,_simulator,_frequency=1):
        SteppableBasePy.__init__(self,_simulator,_frequency)
    def step(self,mcs):
        for cell in self.cellList:
            # if cell.type == self.CANCER1 and np.random.uniform() < 0.0001:
            #     cell.type = self.CANCER2
            cell.targetVolume+=0.05 

                
    # alternatively if you want to make growth a function of chemical concentration uncomment lines below and comment lines above        
        # field=CompuCell.getConcentrationField(self.simulator,"PUT_NAME_OF_CHEMICAL_FIELD_HERE")
        # pt=CompuCell.Point3D()
        # for cell in self.cellList:
            # pt.x=int(cell.xCOM)
            # pt.y=int(cell.yCOM)
            # pt.z=int(cell.zCOM)
            # concentrationAtCOM=field.get(pt)
            # cell.targetVolume+=0.01*concentrationAtCOM  # you can use here any fcn of concentrationAtCOM     
        
        

class MitosisSteppable(MitosisSteppableBase):
    def __init__(self,_simulator,_frequency=1):
        MitosisSteppableBase.__init__(self,_simulator, _frequency)
    
    def step(self,mcs):
        # print "INSIDE MITOSIS STEPPABLE"
        cells_to_divide=[]
        for cell in self.cellList:
            if (cell.type == self.CANCER1 and cell.volume > 30) or cell.volume>35:
                cells_to_divide.append(cell)
                
        for cell in cells_to_divide:
            # to change mitosis mode leave one of the below lines uncommented
            self.divideCellRandomOrientation(cell)
            divide_times['last_division'] = mcs
            # self.divideCellOrientationVectorBased(cell,1,0,0)                 # this is a valid option
            # self.divideCellAlongMajorAxis(cell)                               # this is a valid option
            # self.divideCellAlongMinorAxis(cell)                               # this is a valid option

    def updateAttributes(self):
        parentCell=self.mitosisSteppable.parentCell
        childCell=self.mitosisSteppable.childCell

        divide_times[parentCell.id] = divide_times['last_division']
        divide_times[childCell.id] = divide_times['last_division']
        parentCell.targetVolume = 25 if parentCell.type == 1 else 35
        childCell.targetVolume = 25 if parentCell.type == 1 else 35
        childCell.lambdaVolume = 2
        
        childCell.type = parentCell.type




        

class DeathSteppable(SteppableBasePy):
    def __init__(self,_simulator,_frequency=1):
        SteppableBasePy.__init__(self,_simulator,_frequency)
    def step(self,mcs):
        print divide_times
        for cell in self.cellList:
            if mcs - divide_times[cell.id] > 450:
                cell.targetVolume -= 0.05
                cell.lambdaVolume = 1

        pass
        # if mcs==1000:
        #     for cell in self.cellList:
        #         if cell.type==1:
        #             cell.targetVolume==0
        #             cell.lambdaVolume==100
        
        