
from PySteppables import *
import CompuCell
import sys
import CompuCellSetup
import numpy as np

from PySteppablesExamples import MitosisSteppableBase

phenotype_template = {
    # 'deleterious':(0,0.7*10**15),
    # 15 = genome_order
    'advantageous': (0, 0.1 * 10**15)
}

save_flag = True
simulate_flag = True

divide_times = {'last_division':0}
GLOBAL = {
    'targetVolume':50,
    'divideThreshold':65,
    'cancer2_divideThreshold':64,
    'cancer1_divideThreshold':64,
    'maxTargetVolume':75,
    'cancer2_additional_dV':0.06,
    'cancer1_additional_dV':0.06,
    'init_num_divisions':40,
    'dV':0.05,
    'p_csc_end':0.05,
    'p_csc':0.4
}


import time 
time_info = '_'.join(time.asctime().split(' '))

import os
save_dir = '../simulation_out/CSCs_'+time_info
os.makedirs( save_dir )

from cc3dtools.Tracker import Tracker2
from cc3dtools.Genome import Genome, save_genomes2
from cc3dtools.Phenotype import Phenotype
genomes = {}
phenotypes = {}
divisions_left = {}

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
                # phenotypes[cell.id] = Phenotype( phenotype_template )
                divisions_left[cell.id] = GLOBAL['init_num_divisions']


            if cell.id == r:
                cell.type = self.CANCER2

                if simulate_flag:
                    genomes[cell.id] = Genome( mutation_rate = 120 , name = cell.id,  ploidy_probability=0.002 , ploidy=2 )
                    divisions_left[cell.id] = -1


            if save_flag:
                self.start_tracker.stash( [ cell.id, cell.type , genomes[cell.id].mutation_rate ] )
            # holder[cell.id] = { 'g': Genome( mutation_rate = 20 , genome_order = 10 ) }


        ## check if we successfully initiated a cancer cell.
        count = 0
        for cell in self.cellList:
            if cell.type == self.CANCER2:
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

        y_intercept = np.log( 1/float(GLOBAL['dV']) - 1 )

        for cell in self.cellList:
            # don't grow if your volume is very very large
            # print cell.targetVolume - cell.volume
            # if cell.targetVolume - cell.volume < 1: 
                # print 'too large'
            # else:
            # if cell.type == self.CANCER1 and np.random.uniform() < 0.0001:
            #     cell.type = self.CANCER2
            if cell.type == self.NORMAL:
                cell.targetVolume = min( GLOBAL['dV'] + cell.targetVolume , GLOBAL['maxTargetVolume'] )
            
            if cell.type == self.CANCER1:
                cell.targetVolume += GLOBAL['cancer1_additional_dV']

            if cell.type == self.CANCER2:
                cell.targetVolume += GLOBAL['cancer2_additional_dV']


 
        

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
                if divisions_left[cell.id] != 0:
                    cells_to_divide.append(cell)


                
        for cell in cells_to_divide:
            # to change mitosis mode leave one of the below lines uncommented
            divisions_left[cell.id] -= 1
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
            new_mutations = genomes[parentCell.id].mutate()

            # for mutation in new_mutations:
            #     phenotypes[parentCell.id].evaluate( mutation )

            genomes[childCell.id] = genomes[parentCell.id].replicate( name = childCell.id )
            # phenotypes[childCell.id] = phenotypes[parentCell.id].replicate()

        if parentCell.type == self.NORMAL:
            divisions_left[childCell.id] = divisions_left[parentCell.id]
            childCell.type = parentCell.type
        else:
            # logic for CSCs
            if parentCell.type == self.CANCER2:
                r = np.random.rand()
                if r < GLOBAL['p_csc_end']:
                    # CSC spontaneously looses its capability..
                    parentCell.type = self.CANCER1
                    divisions_left[childCell.id] = 10 # 10 more divisions only.
                    divisions_left[parentCell.id] = 10
                    childCell.type = parentCell.type
                else:
                    childCell.type = self.CANCER1
                    divisions_left[childCell.id] = 10

                if r < GLOBAL['p_csc']:
                    # the child cell will be a CSC!
                    childCell.type = self.CANCER2
                    divisions_left[childCell.id] = -1

            else:
                divisions_left[childCell.id] = 10
                childCell.type = self.CANCER1
            

        




        childCell.lambdaVolume = 1.5
        parentCell.lambdaVolume = 1.5
        # if parentCell.type == self.CANCER1:
        #     childCell.type = self.CANCER2
        # else:
        #     childCell.type = parentCell.type


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

## this demo is to show how to obtain types of the neighbour
class NeighborTrackerPrinterSteppable(SteppableBasePy):
    def __init__(self,_simulator,_frequency=100):
        SteppableBasePy.__init__(self,_simulator,_frequency)

    def step(self,mcs):
        pass
        # for cell in self.cellList:            
            # print "*********NEIGHBORS OF CELL WITH ID ",cell.id," *****************"
        #     count1 = 0
        #     count2 = 0
            # for neighbor , commonSurfaceArea in self.getCellNeighborDataList(cell):                
            #     if neighbor:
            #         if neighbor.type == 2 or neighbor.type == 3:
            #             count2 += 1
            #         else:
            #             count1 += 1

        #             print "neighbor.type",neighbor.type
        #             print "neighbor.id",neighbor.id," commonSurfaceArea=",commonSurfaceArea
        # #         # else:
        #             print "Medium commonSurfaceArea=",commonSurfaceArea
            # print 'non-cancer neighbours:',count1,', cancer neighbours:',count2
        # time.sleep(0.1)


        

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
#         # if mcs==1000:
#         #     for cell in self.cellList:
#         #         if cell.type==1:
#         #             cell.targetVolume==0
#         #             cell.lambdaVolume==100


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
        