
from PySteppables import *
import CompuCell
import sys
import numpy as np
from PySteppablesExamples import MitosisSteppableBase
            

class ConstraintInitializerSteppable(SteppableBasePy):
    def __init__(self,_simulator,_frequency=1):
        SteppableBasePy.__init__(self,_simulator,_frequency)
    def start(self):

        for cell in self.cellList:
            cell.targetVolume=25
            cell.lambdaVolume=2.0
        
        

class GrowthSteppable(SteppableBasePy):
    def __init__(self,_simulator,_frequency=1):
        SteppableBasePy.__init__(self,_simulator,_frequency)
    def step(self,mcs):
        pass
        # for cell in self.cellList:
        #     cell.targetVolume+=1        
        # alternatively if you want to make growth a function of chemical concentration uncomment lines below and comment lines above        
        O2field=CompuCell.getConcentrationField(self.simulator,"OXYGEN")
        MMPfield=CompuCell.getConcentrationField(self.simulator,"MMP")
        pt=CompuCell.Point3D()
        for cell in self.cellList:
            pt.x=int(cell.xCOM)
            pt.y=int(cell.yCOM)
            pt.z=int(cell.zCOM)
            O2Conc=O2field.get(pt)
            MMPConc=MMPfield.get(pt)
            if O2Conc < 0: print('--------~~~~~~~~~-----> O2 VERY LOW') 
            print '------///------>',cell.targetVolume, O2Conc, MMPConc
            if cell.type == self.NORM and MMPConc > 1:
                print 'MMP CONC IS CRAZY!'
                # raw_input('!')
                cell.targetVolume -= 0.5
                cell.lambdaVolume = 2
                # raw_input('!-->removing a NORM cell')
            else:
            # if True:
                if cell.type == self.TPROL and O2Conc < 0:
                    cell.type = self.TMIGR
                    continue
                # O2Conc = np.abs(O2Conc)

                cell.targetVolume+=0.01*O2Conc / (0.05 + O2Conc)  # you can use here any fcn of concentrationAtCOM     
        
        

class MitosisSteppable(MitosisSteppableBase):
    def __init__(self,_simulator,_frequency=1):
        MitosisSteppableBase.__init__(self,_simulator, _frequency)
    
    def step(self,mcs):
        # print "INSIDE MITOSIS STEPPABLE"
        cells_to_divide=[]
        for cell in self.cellList:
            if cell.volume>42 or ( cell.type == self.TPROL and cell.volume > 40 ):
                pass
                # raw_input('celltype:'+str(cell.type))
                cells_to_divide.append(cell)
                
        for cell in cells_to_divide:
            # to change mitosis mode leave one of the below lines uncommented
            self.divideCellRandomOrientation(cell)
            # self.divideCellOrientationVectorBased(cell,1,0,0)                 # this is a valid option
            # self.divideCellAlongMajorAxis(cell)                               # this is a valid option
            # self.divideCellAlongMinorAxis(cell)                               # this is a valid option

    def updateAttributes(self):
        parentCell=self.mitosisSteppable.parentCell
        childCell=self.mitosisSteppable.childCell
        
        childCell.targetVolume=25
        childCell.lambdaVolume=parentCell.lambdaVolume

        
        

class DeathSteppable(SteppableBasePy):
    def __init__(self,_simulator,_frequency=1):
        SteppableBasePy.__init__(self,_simulator,_frequency)
    def step(self,mcs):
        # if mcs==1000:
        #     for cell in self.cellList:
        #         if cell.type==1:
        #             cell.targetVolume==0
        #             cell.lambdaVolume==100
        pass
        


from PySteppables import *
import CompuCellSetup



class PlotSteppable(SteppableBasePy):
    def __init__(self,_simulator,_frequency=10):
        SteppableBasePy.__init__(self,_simulator,_frequency)

    def start(self):
        
        # avg volumes
        self.pWVol=CompuCellSetup.addNewPlotWindow(_title='Average Volume',_xAxisTitle='MonteCarlo Step (MCS)',_yAxisTitle='Average Volume')        
        self.pWVol.addPlot(_plotName='MVol',_style='Dots',_color='red',_size=5)        
        self.pWVol.addPlot(_plotName='N_TVol',_style='Dots',_color='green',_size=5)        


        # number of cells
        self.pWNum=CompuCellSetup.addNewPlotWindow(_title='Number of Cells',_xAxisTitle='MonteCarlo Step (MCS)',_yAxisTitle='Number')                
        self.pWNum.addPlot(_plotName='Normal',_color='green',_size=2)
        self.pWNum.addPlot(_plotName='PROL',_color='blue',_size=2)
        self.pWNum.addPlot(_plotName='MIGR',_color='red', _style='Dots',_size=2)
        self.pWNum.addPlot(_plotName='NECR',_color='black', _style='Dots',_size=2)
        self.pWNum.addPlot(_plotName='Total',_color='red',_size=2)

        # proportions
        self.pWProp=CompuCellSetup.addNewPlotWindow(_title='Proportions of Cancer vs Normal Cells',_xAxisTitle='MonteCarlo Step (MCS)',_yAxisTitle='%')                
        self.pWProp.addPlot(_plotName='Cancer',_color='green')
        self.pWProp.addPlot(_plotName='Normal',_color='blue')


    def step(self,mcs):
        cancer = 0
        necr = 0
        migr = 0
        normal = 0
        meanSurface=0.0
        meanVolume=0.0
        numberOfCells=0
        meanTargetVolume = 0
        for cell  in  self.cellList:
            if cell.type != 0:
                meanVolume+=cell.volume
                meanSurface+=cell.surface
                numberOfCells+=1
                if cell.type == 1:
                    normal+=1
                    meanTargetVolume+=cell.targetVolume
                elif cell.type == 2:
                    cancer+=1
                elif cell.type == 3:
                    migr+=1
                elif cell.type == self.NECR: 
                    necr+=1
                else:
                    pass

        meanVolume/=float(numberOfCells)
        meanTargetVolume/=float(numberOfCells)
        meanSurface/=float(numberOfCells)
        
        self.pWVol.addDataPoint("MVol",mcs,meanVolume)
        self.pWVol.addDataPoint("N_TVol",mcs,meanTargetVolume)

        self.pWNum.addDataPoint("PROL",mcs,cancer)
        self.pWNum.addDataPoint("Normal",mcs,normal)
        self.pWNum.addDataPoint("MIGR",mcs,migr)
        self.pWNum.addDataPoint("NECR",mcs,necr)
        self.pWNum.addDataPoint("Total",mcs,numberOfCells)


        self.pWProp.addDataPoint("Cancer",mcs,cancer/float(numberOfCells))
        self.pWProp.addDataPoint("Normal",mcs,normal/float(numberOfCells))
        

        # print "meanVolume=",meanVolume,"meanSurface=",meanSurface
                
        self.pWVol.showAllPlots()
        self.pWNum.showAllPlots()
        self.pWProp.showAllPlots()
