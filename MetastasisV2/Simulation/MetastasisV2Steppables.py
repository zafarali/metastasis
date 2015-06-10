
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
                cell.targetVolume -= 0.5
                cell.lambdaVolume = 5
                # raw_input('!-->removing a NORM cell')
            else:
                if cell.type == self.TPROL and O2Conc < 0:
                    cell.type = self.TMIGR
                    continue
                O2Conc = np.abs(O2Conc)
                cell.targetVolume+=0.01*O2Conc / (0.05 + O2Conc)  # you can use here any fcn of concentrationAtCOM     
        
        

class MitosisSteppable(MitosisSteppableBase):
    def __init__(self,_simulator,_frequency=1):
        MitosisSteppableBase.__init__(self,_simulator, _frequency)
    
    def step(self,mcs):
        # print "INSIDE MITOSIS STEPPABLE"
        cells_to_divide=[]
        for cell in self.cellList:
            if cell.volume>42:
                
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
        if mcs==1000:
            for cell in self.cellList:
                if cell.type==1:
                    cell.targetVolume==0
                    cell.lambdaVolume==100
        
        