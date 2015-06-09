

from PySteppables import *
import CompuCellSetup



class PlotSteppable(SteppableBasePy):
    def __init__(self,_simulator,_frequency=10):
        SteppableBasePy.__init__(self,_simulator,_frequency)

    def start(self):
        
        # avg volumes
        self.pWVol=CompuCellSetup.addNewPlotWindow(_title='Average Volume',_xAxisTitle='MonteCarlo Step (MCS)',_yAxisTitle='Average Volume')        
        self.pWVol.addPlot(_plotName='MVol',_style='Dots',_color='red',_size=5)        


        # number of cells
        self.pWNum=CompuCellSetup.addNewPlotWindow(_title='Number of Cells',_xAxisTitle='MonteCarlo Step (MCS)',_yAxisTitle='Number')                
        self.pWNum.addPlot(_plotName='Cancer',_color='green')
        self.pWNum.addPlot(_plotName='Normal',_color='blue')
        self.pWNum.addPlot(_plotName='Total',_color='red')

        # proportions
        self.pWProp=CompuCellSetup.addNewPlotWindow(_title='Proportions of Cancer vs Normal Cells',_xAxisTitle='MonteCarlo Step (MCS)',_yAxisTitle='%')                
        self.pWProp.addPlot(_plotName='Cancer',_color='green')
        self.pWProp.addPlot(_plotName='Normal',_color='blue')


    def step(self,mcs):
        cancer = 0
        normal = 0
        meanSurface=0.0
        meanVolume=0.0
        numberOfCells=0
        for cell  in  self.cellList:
            meanVolume+=cell.volume
            meanSurface+=cell.surface
            numberOfCells+=1
            if cell.type == 1:
                normal+=1
            else:
                cancer+=1

        meanVolume/=float(numberOfCells)
        meanSurface/=float(numberOfCells)
        
        self.pWVol.addDataPoint("MVol",mcs,meanVolume)
        self.pWNum.addDataPoint("Cancer",mcs,cancer)
        self.pWNum.addDataPoint("Normal",mcs,normal)
        self.pWNum.addDataPoint("Total",mcs,numberOfCells)


        self.pWProp.addDataPoint("Cancer",mcs,cancer/float(numberOfCells))
        self.pWProp.addDataPoint("Normal",mcs,normal/float(numberOfCells))
        

        # print "meanVolume=",meanVolume,"meanSurface=",meanSurface
                
        self.pWVol.showAllPlots()
        self.pWNum.showAllPlots()
        self.pWProp.showAllPlots()
