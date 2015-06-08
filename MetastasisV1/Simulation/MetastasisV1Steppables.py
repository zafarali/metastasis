
from PySteppables import *
import CompuCell
import sys

import numpy as np

from PySteppablesExamples import MitosisSteppableBase
            
# phenotypes = {
#     'deleterious':(0,0.7)
# }

divide_times = {'last_division':0}

from cc3dtools.Tracker import Tracker2
from cc3dtools.Genome import Genome, save_genomes
genomes = {}

class ConstraintInitializerSteppable(SteppableBasePy):
    def __init__(self,_simulator,_frequency=1):
        SteppableBasePy.__init__(self,_simulator,_frequency)
    def start(self):
        start_tracker = Tracker2(file_name='../start_cells_862015.csv')

        num_cells = len(self.cellList)

        r = np.random.randint(0,num_cells)

        for cell in self.cellList:
            
            divide_times[cell.id] = 0
            
            cell.targetVolume=50
            cell.lambdaVolume=1

            genomes[cell.id] = Genome( mutation_rate = 50 , name = cell.id )

            if cell.id == r:
                cell.type = self.CANCER1
                genomes[cell.id] = Genome( mutation_rate = 120 , name = cell.id )

            start_tracker.stash( [ cell.id, cell.type , genomes[cell.id].mutation_rate ] )
            # holder[cell.id] = { 'g': Genome( mutation_rate = 20 , genome_order = 10 ), 'p': Phenotype() }

        start_tracker.save_stash(flag='w')

    def finish(self):
        
        tracker = Tracker2(file_name='../finish_cells_862015.csv')

        for cell in self.cellList:
            z = cell.zCOM
            y = cell.yCOM
            x = cell.xCOM
            tracker.stash( [ cell.id , cell.type , x , y , z ] )

        # tracker.save_stash()
        # save_genomes( [ genome[1] for genome in genomes.items() ] , file_name = '../genomes_862015.csv', method = 'aligned' )

        

        

class GrowthSteppable(SteppableBasePy):
    def __init__(self,_simulator,_frequency=1):
        SteppableBasePy.__init__(self,_simulator,_frequency)
    def step(self,mcs):
        for cell in self.cellList:
            # don't grow if your volume is very very large
            # if cell.targetVolume - cell.volume > 10: continue
            # if cell.type == self.CANCER1 and np.random.uniform() < 0.0001:
            #     cell.type = self.CANCER2
            cell.targetVolume+=0.05 
            if cell.type == self.CANCER1 or cell.type == self.CANCER2:
                # cancerous cells grow slightly faster
                cell.targetVolume += 0.05
            # if cell.type == self.CANCER2:
                # cell.targetVolume += 1
                # cell.lambdaVolume = 3

                
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
    def __init__(self,_simulator, tracker_instance ,_frequency=1):
        MitosisSteppableBase.__init__(self,_simulator, _frequency)
        self.tracker_instance = tracker_instance
    def start(self):
        # we initialize the stash function
        for cell in self.cellList:
                                            # R for 'root'
            self.tracker_instance.stashDivision( 'R' , cell.id, cell.id, cell.id )


    def step(self,mcs):
        # print "INSIDE MITOSIS STEPPABLE"
        cells_to_divide=[]
        for cell in self.cellList:
            if (cell.type == self.CANCER1 and cell.volume > 80) or cell.volume>100:
            # if cell.volume > 100:
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
        childCell.targetVolume=50
        parentCell.targetVolume=50

        # if parentCell.type == self.CANCER1:
        #     childCell.type = self.CANCER1
        #     childCell.targetVolume = 40
        genomes[parentCell.id].mutate()
        genomes[childCell.id] = genomes[parentCell.id].replicate( name = childCell.id )

        childCell.lambdaVolume = 1.5
        parentCell.lambdaVolume = 1.5
        if parentCell.type == self.CANCER1:
            childCell.type = self.CANCER2
        else:
            childCell.type = parentCell.type

        self.tracker_instance.stashDivision( divide_times['last_division'] , parentCell.id, childCell.id, parentCell.id )

    def finish(self):
        self.tracker_instance.saveStash()
        

class DeathSteppable(SteppableBasePy):
    def __init__(self,_simulator,_frequency=1):
        SteppableBasePy.__init__(self,_simulator,_frequency)
    def step(self,mcs):
        # print divide_times
        # for cell in self.cellList:
        #     if mcs - divide_times[cell.id] > 450:
        #         cell.targetVolume = 0 
        #         cell.lambdaVolume = 100

        pass
#         # if mcs==1000:
#         #     for cell in self.cellList:
#         #         if cell.type==1:
#         #             cell.targetVolume==0
#         #             cell.lambdaVolume==100
        
