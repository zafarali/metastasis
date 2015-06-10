
from PySteppables import *
import CompuCell
import sys

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
        field=CompuCell.getConcentrationField(self.simulator,"OXYGEN")
        pt=CompuCell.Point3D()
        for cell in self.cellList:
            pt.x=int(cell.xCOM)
            pt.y=int(cell.yCOM)
            pt.z=int(cell.zCOM)
            concentrationAtCOM=field.get(pt)
            if concentrationAtCOM < 0: raw_input('!')
            print '------///------>',cell.targetVolume, concentrationAtCOM
            cell.targetVolume+=0.01*concentrationAtCOM / (0.05 + concentrationAtCOM)  # you can use here any fcn of concentrationAtCOM     
        
        

class MitosisSteppable(MitosisSteppableBase):
    def __init__(self,_simulator,_frequency=1):
        MitosisSteppableBase.__init__(self,_simulator, _frequency)
    
    def step(self,mcs):
        # print "INSIDE MITOSIS STEPPABLE"
        cells_to_divide=[]
        for cell in self.cellList:
            if cell.volume>50:
                
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
        
        childCell.targetVolume=parentCell.targetVolume
        childCell.lambdaVolume=parentCell.lambdaVolume
        if parentCell.type==1:
            childCell.type=2
        else:
            childCell.type=1
        
        

class DeathSteppable(SteppableBasePy):
    def __init__(self,_simulator,_frequency=1):
        SteppableBasePy.__init__(self,_simulator,_frequency)
    def step(self,mcs):
        if mcs==1000:
            for cell in self.cellList:
                if cell.type==1:
                    cell.targetVolume==0
                    cell.lambdaVolume==100
        
        