
from PySteppables import *
import CompuCell
import sys
import CompuCellSetup
import numpy as np

from PySteppablesExamples import MitosisSteppableBase


save_flag = True
simulate_flag = True

divide_times = {'last_division':0}
GLOBAL = {
    'targetVolume':50,
    'divideThreshold':65,
    'cancer2_divideThreshold':63,
    'cancer1_divideThreshold':63,
    'maxTargetVolume':75,
    'cancer2_additional_dV':0.03,
    'cancer1_additional_dV':0.03,
    'dV':0.05
}

import time 
time_info = '_'.join(time.asctime().split(' '))

import os
save_dir = '../simulation_out/minor_adv_'+time_info
os.makedirs( save_dir )


from cc3dtools.Tracker import Tracker2
from cc3dtools.Genome import Genome, save_genomes2
genomes = {}

class ConstraintInitializerSteppable(SteppableBasePy):
    def __init__(self,_simulator,_frequency=1):
        SteppableBasePy.__init__(self,_simulator,_frequency)
        if save_flag:
            self.start_tracker = Tracker2( file_name = save_dir+'/start_cells_'+time_info+'.csv')

    def start(self):
        
        num_cells = len(self.cellList)

        r = np.random.randint(0,num_cells)

        for cell in self.cellList:
            
            divide_times[cell.id] = 0
            
            cell.targetVolume=GLOBAL['targetVolume']
            cell.lambdaVolume=1.5

            if simulate_flag:
                genomes[cell.id] = Genome( mutation_rate = 50 , name = cell.id, ploidy_probability=0.0001 , ploidy=2 )

            if cell.id == r:
                cell.type = self.CANCER1

                if simulate_flag:
                    genomes[cell.id] = Genome( mutation_rate = 120 , name = cell.id, ploidy_probability=0.002 , ploidy=2 )

            if save_flag:
                self.start_tracker.stash( [ cell.id, cell.type , genomes[cell.id].mutation_rate ] )
            # holder[cell.id] = { 'g': Genome( mutation_rate = 20 , genome_order = 10 ) }


        ## check if we successfully initiated a cancer cell.
        count = 0
        for cell in self.cellList:
            if cell.type == self.CANCER1:
                count += 1

        if count != 1:
            print 'FAILED TO INITIALZE A CANCER CELL.'
            sys.exit(0)

        

    def finish(self):
        if save_flag:
            tracker = Tracker2( file_name = save_dir+'/finish_cells_'+time_info+'.csv' )

            for cell in self.cellList:
                z = cell.zCOM
                y = cell.yCOM
                x = cell.xCOM
                tracker.stash( [ cell.id , cell.type , x , y , z ] )

            tracker.save_stash() # save final cell data
            self.start_tracker.save_stash() # save initial cell data
            save_genomes2( [ genome[1] for genome in genomes.items() ] , file_name = save_dir+'/genomes_'+time_info+'.csv' ) #save genomes

        

        

class GrowthSteppable(SteppableBasePy):
    def __init__(self,_simulator,_frequency=1):
        SteppableBasePy.__init__(self,_simulator,_frequency)
    def step(self,mcs):
        #                   print '----->Global targetVolume',GLOBAL['targetVolume']
        for cell in self.cellList:
            # don't grow if your volume is very very large
            # print cell.targetVolume - cell.volume
            # if cell.targetVolume - cell.volume < 1: 
                # print 'too large'
            # else:
            # if cell.type == self.CANCER1 and np.random.uniform() < 0.0001:
            #     cell.type = self.CANCER2
            cell.targetVolume = min( GLOBAL['dV'] + cell.targetVolume , GLOBAL['maxTargetVolume'] )
            
            # cell.targetVolume = min(0.05 + cell.targetVolume if cell.targetVolume < GLOBAL['divideThreshold'] + 3 else cell.targetVolume, GLOBAL['maxTargetVolume'])

            #                       print '------> growth event:',cell.targetVolume, cell.volume
            # if cell.type == self.CANCER1 or cell.type == self.CANCER2:
            if cell.type == self.CANCER2:
                # cancerous cells grow slightly faster
                cell.targetVolume += GLOBAL['cancer2_additional_dV']

                #               print '------>cancer growth event:',cell.targetVolume, cell.volume

            if cell.type == self.CANCER1:
                # cancerous cells grow slightly faster
                cell.targetVolume += GLOBAL['cancer1_additional_dV']
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
    def __init__(self,_simulator,_frequency=1):
        MitosisSteppableBase.__init__(self,_simulator, _frequency)
        if save_flag:
            self.mitosis_tracker = Tracker2(file_name=save_dir+'/division_events_'+time_info+'.csv')

    def start(self):
        # we initialize the stash function
        if save_flag:
            for cell in self.cellList:
                                        # R for 'root'
                self.mitosis_tracker.stash( [ 'R' , cell.id, cell.id, cell.id ] )


    def step(self,mcs):
        # print "INSIDE MITOSIS STEPPABLE"
        cells_to_divide=[]
        for cell in self.cellList:
            if ( cell.type == self.CANCER2 and cell.volume > GLOBAL['cancer2_divideThreshold'] ) or \
            ( cell.type == self.CANCER1 and cell.volume > GLOBAL['cancer1_divideThreshold'] ) or \
            cell.volume > GLOBAL['divideThreshold'] :
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
        #                           print '------>mitosis event:',childCell.volume, parentCell.volume
        # raw_input()
        divide_times[parentCell.id] = divide_times['last_division']
        divide_times[childCell.id] = divide_times['last_division']

        ## using the equation derived to keep force in the model constant
        ## T_star = T_target - t/2 <-- t/2 is volume after division 
        T_star = GLOBAL['targetVolume'] - ( parentCell.volume / 2.0 )


        childCell.targetVolume = T_star
        parentCell.targetVolume = T_star
        # GLOBAL['targetVolume'] = GLOBAL['targetVolume'] - 0.05 if GLOBAL['targetVolume'] > 45 else GLOBAL['targetVolume']
        # if parentCell.type == self.CANCER1:
        #     childCell.type = self.CANCER1
        #     childCell.targetVolume = 40

        if simulate_flag:
            genomes[parentCell.id].mutate()
            genomes[childCell.id] = genomes[parentCell.id].replicate( name = childCell.id )

        childCell.lambdaVolume = 1.5
        parentCell.lambdaVolume = 1.5
        # if parentCell.type == self.CANCER1:
        #     childCell.type = self.CANCER2
        # else:
        #     childCell.type = parentCell.type

        childCell.type = parentCell.type

        # # attempt to obtain the proliferating front
        # if parentCell.type != self.NORMAL:
        #     for cell in [ childCell , parentCell ]:
        #         normal_count = 0
        #         for neighbor , commonSurfaceArea in self.getCellNeighborDataList(cell):                
        #             if neighbor:
        #                 if neighbor.type == self.NORMAL:
        #                     normal_count += 1
        #                 #endif
        #             #endif
        #         cell.type = self.CANCER2 if normal_count > 0 else self.CANCER1
        #     #endfor
        # else:
        #     childCell.type = parentCell.type
        # #endelse

        if save_flag:
            self.mitosis_tracker.stash( [ divide_times['last_division'] , parentCell.id, childCell.id, parentCell.id ] )

    def finish(self):
        if save_flag:
            self.mitosis_tracker.save_stash()
        pass

class DeathSteppable(SteppableBasePy):
    def __init__(self,_simulator,_frequency=1):
        SteppableBasePy.__init__(self,_simulator,_frequency)

    def step(self,mcs):
        # print divide_times
        # for cell in self.cellList:
        #     if mcs - divide_times[cell.id] > 450:
        #         cell.targetVolume -= 0.007 
        #         cell.lambdaVolume = 1

        pass

class SuperTracker(SteppableBasePy):
    def __init__(self,_simulator,_frequency=10):
        SteppableBasePy.__init__(self,_simulator,_frequency)
        self.cell_tracker = Tracker2( file_name = save_dir+'/cell_count_'+time_info+'.csv' )
        self.volume_tracker = Tracker2( file_name = save_dir+'/volume_'+time_info+'.csv' )

    def step(self,mcs):
        cancer1 = 0
        cancer2 = 0
        normal = 0
        mean_normal_volume = 0
        mean_cancer1_volume = 0
        mean_cancer2_volume = 0
        all_cells = 0
        mean_overall_volume = 0

        for cell in self.cellList:
            if cell.type == self.NORMAL:
                mean_normal_volume += cell.volume
                normal += 1
            elif cell.type == self.CANCER1:
                mean_cancer1_volume += cell.volume
                cancer1 += 1
            elif cell.type == self.CANCER2:
                mean_cancer2_volume += cell.volume
                cancer2 += 1
            
            all_cells += 1
            mean_overall_volume += cell.volume

        mean_overall_volume = mean_overall_volume/float(all_cells) if all_cells > 0 else 0
        mean_normal_volume = mean_normal_volume/float(normal) if normal > 0 else 0
        mean_cancer2_volume =mean_cancer2_volume/ float(cancer2) if cancer2 > 0 else 0
        mean_cancer1_volume =mean_cancer1_volume/ float(cancer1) if cancer1 > 0 else 0

        self.cell_tracker.stash( [ mcs , 0 , normal ] )
        self.cell_tracker.stash( [ mcs , 1 , cancer1 ] )
        self.cell_tracker.stash( [ mcs , 2 , cancer2 ] )
        self.cell_tracker.stash( [ mcs , -1 , all_cells ] )


        self.volume_tracker.stash( [ mcs , 0 , mean_normal_volume ] )
        self.volume_tracker.stash( [ mcs , 1 , mean_cancer1_volume ] )
        self.volume_tracker.stash( [ mcs , 2 , mean_cancer2_volume ] )
        self.volume_tracker.stash( [ mcs , -1 , mean_overall_volume ] )

        
    def finish(self):
        self.cell_tracker.save_stash()
        self.volume_tracker.save_stash()
        