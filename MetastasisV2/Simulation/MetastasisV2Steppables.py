from PySteppables import *
import CompuCell
import sys

from PySteppablesExamples import MitosisSteppableBase

import numpy as np            

class ConstraintInitializerSteppable(SteppableBasePy):
    def __init__(self,_simulator,_frequency=1):
        SteppableBasePy.__init__(self,_simulator,_frequency)
    def start(self):
        for cell in self.cellList:
            cell.targetVolume=35
            cell.lambdaVolume=2.0
        
        

class GrowthSteppable(SteppableBasePy):
    def __init__(self,_simulator,_frequency=1):
        SteppableBasePy.__init__(self,_simulator,_frequency)
    def step(self,mcs):
        for cell in self.cellList:
            cell.targetVolume+=0.1
        #     if cell.type == self.TPROL:
        #         cell.targetVolume += 0.05
        # #     # pass
        # alternatively if you want to make growth a function of chemical concentration uncomment lines below and comment lines above        
        # field=CompuCell.getConcentrationField(self.simulator,"OXYGEN")
        # pt=CompuCell.Point3D()
        # for cell in self.cellList:
        #     if cell.type == self.WALL: continue
        #     pt.x=int(cell.xCOM)
        #     pt.y=int(cell.yCOM)
        #     pt.z=int(cell.zCOM)
        #     concentrationAtCOM= np.abs( field.get(pt) )
        #     cell.targetVolume+=0.001*concentrationAtCOM  # you can use here any fcn of concentrationAtCOM     
        
        

class MitosisSteppable(MitosisSteppableBase):
    def __init__(self,_simulator,_frequency=1):
        MitosisSteppableBase.__init__(self,_simulator, _frequency)
    
    def step(self,mcs):
        # print "INSIDE MITOSIS STEPPABLE"
        cells_to_divide=[]
        for cell in self.cellList:
            if cell.volume>80 and cell.type != self.WALL:
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
        parentCell.targetVolume=35
        childCell.targetVolume=35
        childCell.lambdaVolume=parentCell.lambdaVolume
        # if parentCell.type==1:
        #     childCell.type=2
        # else:
        #     childCell.type=1
        
        

class DeathSteppable(SteppableBasePy):
    def __init__(self,_simulator,_frequency=100):
        SteppableBasePy.__init__(self,_simulator,_frequency)
    def step(self,mcs):
        for cell in self.cellList:

            OxField = CompuCell.getConcentrationField( self.simulator , "OXYGEN" )
            MMPField = CompuCell.getConcentrationField( self.simulator , "MMP" )
            pt = CompuCell.Point3D()
            pt.x = int(cell.xCOM)
            pt.y = int(cell.yCOM)
            pt.z = int(cell.zCOM)

            O2conc = OxField.get( pt )
            print '--->O2',O2conc

            MMPconc = MMPField.get( pt )
            print '--->MMP',MMPconc

            if MMPconc < 0 or O2conc < 0:
                raw_input('!')
            # if MMPconc < 1e-1 and cell.type == self.NORM:
            #     cell.targetVolume=0
            #     cell.lambdaVolume=100
            if O2conc < 0 and cell.type == self.TPROL:

                cell.type = self.TMIGR



                # if O2conc < 
        # for cell in 
        # if mcs==1000:
        #     for cell in self.cellList:
        #         if cell.type==1:
        #             cell.targetVolume==0
        #             cell.lambdaVolume==100
        pass
        
        