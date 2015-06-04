
from PySteppables import *
import CompuCell
import sys

import numpy as np

from PySteppablesExamples import MitosisSteppableBase
            
# phenotypes = {
#     'deleterious':(0,0.7)
# }

divide_times = {'last_division':0}

class ConstraintInitializerSteppable(SteppableBasePy):
    def __init__(self,_simulator,_frequency=1):
        SteppableBasePy.__init__(self,_simulator,_frequency)
    def start(self):

        num_cells = len(self.cellList)


        for cell in self.cellList:
            
            divide_times[cell.id] = 0

            cell.targetVolume=50
            cell.lambdaVolume=1
            # holder[cell.id] = { 'g': Genome( mutation_rate = 20 , genome_order = 10 ), 'p': Phenotype() }

        r = np.random.randint(0,num_cells)
        self.cellList.inventory.attemptFetchingCellById(r).type = self.CANCER1
        self.cellList.inventory.attemptFetchingCellById(r).targetVolume = 50
        

        

class GrowthSteppable(SteppableBasePy):
    def __init__(self,_simulator,_frequency=1):
        SteppableBasePy.__init__(self,_simulator,_frequency)
    def step(self,mcs):
        for cell in self.cellList:
            # don't grow if your volume is very very large
            if cell.targetVolume - cell.volume > 10: continue
            # if cell.type == self.CANCER1 and np.random.uniform() < 0.0001:
            #     cell.type = self.CANCER2
            cell.targetVolume+=0.05 
            if cell.type == self.CANCER1:
                # cancerous cells grow slightly faster
                cell.targetVolume += 0.025

                
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
            # if (cell.type == self.CANCER1 and cell.volume > 30) or cell.volume>35:
            if cell.volume > 100:
                cells_to_divide.append(cell)
                
        for cell in cells_to_divide:
            # to change mitosis mode leave one of the below lines uncommented
            # self.divideCellRandomOrientation(cell)
            divide_times['last_division'] = mcs
            # self.divideCellOrientationVectorBased(cell,1,0,0)                 # this is a valid option
            # self.divideCellAlongMajorAxis(cell)                               # this is a valid option
            self.divideCellAlongMinorAxis(cell)                               # this is a valid option

    def updateAttributes(self):
        parentCell=self.mitosisSteppable.parentCell
        childCell=self.mitosisSteppable.childCell

        divide_times[parentCell.id] = divide_times['last_division']
        divide_times[childCell.id] = divide_times['last_division']
        childCell.targetVolume=50
        childCell.lambdaVolume=50
        # if parentCell.type == self.CANCER1:
        #     childCell.type = self.CANCER1
        #     childCell.targetVolume = 40

        childCell.lambdaVolume = 1
        
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
        
      
# class ExtraMultiPlotSteppable(SteppablePy):
#     def __init__(self,_simulator,_frequency=10):
#         SteppablePy.__init__(self,_frequency)
#         self.simulator=_simulator
#         self.inventory=self.simulator.getPotts().getCellInventory()
#         self.cellList=CellList(self.inventory)

#     def start(self):
#         # import CompuCellSetup  
#         # print "CompuCellSetup.viewManager=",CompuCellSetup.viewManager    
#         # CompuCellSetup.viewManager.plotManager.addNewPlotWindow()

#         import CompuCellSetup  
#         #self.pW=CompuCellSetup.addNewPlotWindow()
#         # CompuCellSetup.viewManager.plotManager.emitNewPlotWindow()
        
#         self.pWVol=CompuCellSetup.viewManager.plotManager.getNewPlotWindow()    
#         if not self.pWVol:
#             return         
#         self.pWVol.setTitle("Average Volume")        
#         self.pWVol.setXAxisTitle("MonteCarlo Step (MCS)")
#         self.pWVol.setYAxisTitle("Average Volume")        
#         self.pWVol.addPlot("MVol")
#         self.pWVol.changePlotProperty("MVol","LineWidth",5)
#         self.pWVol.changePlotProperty("MVol","LineColor","red")     
#         self.pWVol.addGrid()
#         #adding automatically generated legend
#         self.pWVol.addAutoLegend()
        
#         self.pWSur=CompuCellSetup.viewManager.plotManager.getNewPlotWindow()        
#         self.pWSur.setTitle("Average Surface")        
#         self.pWSur.setXAxisTitle("MonteCarlo Step (MCS)")        
#         self.pWSur.setYAxisTitle("Average Surface")        
#         self.pWSur.addPlot("MSur")
#         self.pWSur.changePlotProperty("MSur","LineWidth",1)
#         self.pWSur.changePlotProperty("MSur","LineColor","green")         
#         self.pWSur.addGrid()
        
#     def step(self,mcs):
#         # this is totally non optimized code. It is for illustrative purposes only. 
#         meanSurface=0.0
#         meanVolume=0.0
#         numberOfCells=0
#         for cell  in  self.cellList:
#             meanVolume+=cell.volume
#             meanSurface+=cell.surface
#             numberOfCells+=1
#         meanVolume/=float(numberOfCells)
#         meanSurface/=float(numberOfCells)
        
#         self.pWVol.addDataPoint("MVol",mcs,meanVolume)
#         self.pWSur.addDataPoint("MSur",mcs,meanSurface)
#         print "meanVolume=",meanVolume,"meanSurface=",meanSurface
        
#         # self.pW.showPlot("MCS1")
#         self.pWVol.showAllPlots()
#         self.pWSur.showAllPlots()        