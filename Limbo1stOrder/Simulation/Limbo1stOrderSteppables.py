
from PySteppables import *
import CompuCell
import sys
import CompuCellSetup
import numpy as np

from PySteppablesExamples import MitosisSteppableBase

## FIXED PATHS
sys.path.append('/home/zahmed/summer15/metastasis/')

save_flag = True
simulate_flag = True
template_flag = True

TEMPLATE_ROOT = '/home/zahmed/summer15/metastasis/LIMBO_TEMPLATE' # in the future get this from a config file.
TEMPLATES = {
    'start_tracker': TEMPLATE_ROOT + '/start_cells.csv',
    'mitosis_tracker': TEMPLATE_ROOT + '/division_events.csv',
    'cell_tracker': TEMPLATE_ROOT + '/cell_count.csv',
    'volume_tracker': TEMPLATE_ROOT + '/volume.csv',
    'genome_file': TEMPLATE_ROOT + '/genomes.csv.gen2'
}


divide_times = {'last_division':0}
GLOBAL = {
    'targetVolume':50,
    'divideThreshold':45,
    'cancer2_divideThreshold':42,
    'cancer1_divideThreshold':42,
    'maxTargetVolume':75,
    'cancer2_additional_dV':0.1,
    'cancer1_additional_dV':0.1,
    'dV':0, # 0 because this is after the normal cell have grown completely, we do not need them to grow again
    '_dV':0.2
}

LATTICE = {
    'x_max':900,
    'y_max':900,
    'x_min':100,
    'y_min':100,
    'center_x_max':555,
    'center_x_min':545,
    'center_y_max':555,
    'center_y_min':545
}

import time 
time_info = '_'.join(time.asctime().split(' '))

with open(GLOBAL_PATH+'outs.txt', 'a') as f:
    f.write(time_info)
    f.write('\n')

import os
save_dir = '/home/zahmed/summer15/metastasis/simulation_out/limbo_1st'

from cc3dtools.Tracker import Tracker2
from cc3dtools.Genome import Genome, save_genomes2, load_genomes_into_dict

genomes = {}

class ConstraintInitializerSteppable(SteppableBasePy):
    def __init__(self,_simulator,_frequency=1):
        SteppableBasePy.__init__(self,_simulator,_frequency)
        if save_flag:
            if template_flag:
                self.start_tracker = Tracker2( file_name = save_dir+'/start_cells_'+time_info+'.csv' , template = TEMPLATES['start_tracker'] )
                print 'loaded start template'
            else:
                self.start_tracker = Tracker2( file_name = save_dir+'/start_cells_'+time_info+'.csv')

    def start(self):
        
        num_cells = len(self.cellList)
        
        cancer_cell_created = False

        if simulate_flag and template_flag:
            print 'loading genomes'
            genomes.update( load_genomes_into_dict( file_name = TEMPLATES['genome_file'] ) )
            pass

        for cell in self.cellList:
            
            divide_times[cell.id] = 0
            
            cell.targetVolume=GLOBAL['targetVolume']
            cell.lambdaVolume=1.5

            if simulate_flag and not template_flag:
                genomes[cell.id] = Genome( mutation_rate = 50 , name = cell.id, ploidy_probability=0.0001 , ploidy=2 )

            # create a cancer cell in the middle
            if not cancer_cell_created and ( cell.xCOM <= LATTICE['center_x_max'] and cell.xCOM >= LATTICE['center_x_min'] ) and ( cell.yCOM <= LATTICE['center_y_max'] and cell.yCOM >= LATTICE['center_y_min'] ):
                cell.type = self.CANCER1
                cancer_cell_created = True

            if cell.type == self.CANCER1:

                if simulate_flag and not template_flag:
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

        

class UtillitySteppable(SteppableBasePy):
    def __init__(self, _simulator, _frequency=1):
        SteppableBasePy.__init__(self,_simulator,_frequency)

    def step(self, mcs):
        print 'INSIDE UtillitySteppable'
        sys.stdout.flush() # flush the buffer, this allows continous saving into a file
        # if mcs == 10:            
            # self.changeNumberOfWorkNodes(2)
            # print "NUMBER OF WORK NODES INCREASED"


    def finish(self):
        print 'THIS IS OVEA!'



class GrowthSteppable(SteppableBasePy):
    def __init__(self,_simulator,_frequency=1):
        SteppableBasePy.__init__(self,_simulator,_frequency)

    def step(self,mcs):
        print "INSIDE GROWTH STEPPABLE"
        #                   print '----->Global targetVolume',GLOBAL['targetVolume']

        for cell in self.cellList:
            z = cell.zCOM
            y = cell.yCOM
            x = cell.xCOM


            if cell.type == self.DEAD:
                # do not grow if dead
                continue
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
            if template_flag:
                self.mitosis_tracker = Tracker2( file_name=save_dir+'/division_events_'+time_info+'.csv' , template=TEMPLATES['mitosis_tracker'] )
                print 'loaded mitosis_tracker'
            else:
                self.mitosis_tracker = Tracker2( file_name=save_dir+'/division_events_'+time_info+'.csv' )

    def start(self):
        # we initialize the stash function
        if save_flag:
            for cell in self.cellList:
                                        # R for 'root'
                self.mitosis_tracker.stash( [ 'R' , cell.id, cell.id, cell.id ] )


    def step(self,mcs):
        print "INSIDE MITOSIS STEPPABLE"
        cells_to_divide=[]
        for cell in self.cellList:

            if cell.type == self.DEAD:
                # do not divide if dead
                continue


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
        print "---->DIVISION EVENT"
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


class SuperTracker(SteppableBasePy):
    def __init__(self,_simulator,_frequency=100):
        SteppableBasePy.__init__(self,_simulator,_frequency)
        if save_flag:
            if template_flag:
                self.cell_tracker = Tracker2( file_name = save_dir+'/cell_count_'+time_info+'.csv' , template = TEMPLATES['cell_tracker'] )
                self.volume_tracker = Tracker2( file_name = save_dir+'/volume_'+time_info+'.csv' , template = TEMPLATES['volume_tracker'] )
                print 'cell and volume trackers'
            else:
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

        if all_cells==0:
            sys.exit('(!) EXIT DUE TO NO CELLS REMAINING')
        # if cancer1 == 0 and cancer2 == 0:
        #     sys.exit('(!) ALL CANCER CELLS DIED')
        mean_overall_volume = mean_overall_volume/float(all_cells) if all_cells > 0 else 0
        mean_normal_volume = mean_normal_volume/float(normal) if normal > 0 else 0
        mean_cancer2_volume =mean_cancer2_volume/ float(cancer2) if cancer2 > 0 else 0
        mean_cancer1_volume =mean_cancer1_volume/ float(cancer1) if cancer1 > 0 else 0

        print 'MCS@',mcs,' NUMBER OF CELLS REPORT:'
        print 'MCS@',mcs,' Normal cells:',normal
        print 'MCS@',mcs,' cancer1 cells',cancer1
        print 'MCS@',mcs,' cancer2 cells',cancer2
        print 'MCS@',mcs,' all_cells',all_cells

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
        


class DeathCheckSteppable(SteppableBasePy):
    def __init__(self,_simulator,_frequency=10):
        SteppableBasePy.__init__(self,_simulator,_frequency)

    def step(self,mcs):
        print "INSIDE DEATH CHECK STEPPABLE"
        for cell in self.cellList:
            z = cell.zCOM
            y = cell.yCOM
            x = cell.xCOM
            if cell.type == self.DEAD: continue

            if not ( ( x >= LATTICE['x_min'] and x <= LATTICE['x_max'] ) and ( y >= LATTICE['y_min'] and y <=LATTICE['y_max'] ) ) :
                if cell.type == self.CANCER1 or cell.type == self.CANCER2:
                    # self.stopSimulation()
                    pass
                # cell.lambdaVolume = 5
                cell.type = self.DEAD
                self.lambdaVolume = 1
                del genomes[cell.id]
                # print('-------------------------------->LIMBO ZONE')
                # print cell.targetVolume, cell.volume
            
class DeathSteppable(SteppableBasePy):
    def __init__(self, _simulator, _frequency=1):
        SteppableBasePy.__init__(self,_simulator,_frequency)

    def step(self,mcs):
        print "INSIDE DEATH STEPPABLE"
        for cell in self.cellList:
            if cell.type == self.DEAD:
                cell.targetVolume -= 10 * GLOBAL['_dV']

