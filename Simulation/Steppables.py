##########################################################
#	File: MainProgram.py
#	Author: Jose Perez <josegperez@mail.com>
#	Version: Model v5
#
#	This is where we define our steppables.
#	Please see the CC3D documentation for more information about steppables.
#	
##########################################################
from PySteppables import *
from PySteppablesExamples import MitosisSteppableBase
import CompuCell
import sys
import random

# --== Project imports ==--
from Cell import CC3DKey
from Cell import Data
import Config

##########################################################
#	MainSteppable
#
#	This steppable sets-up of our code.
#	We are basically attaching information to the CC3D simulation.
#	That way we can modify and specify our own interactions between cells. 
#
##########################################################
class MainSteppable(SteppableBasePy):   

	### Constructor of our MainSteppable
	###	More info about constructors can be found in the Python documentation
    def __init__(self, _simulator, _frequency = 1):
		# Call our parent's constructor
        SteppableBasePy.__init__(self, _simulator, _frequency)
		
		# Our code doesn't have access to the properties inside this class
		# However we want to be able to change cell types in our code
		# This is sort of a hack that allows that
        # Set-up constants to match CC3D's internal constants
        from Cell import CC3DType
        
        CC3DType.APC = self.APC
        
        CC3DType.TREG_INACTIVE = self.TREG_INACTIVE
        CC3DType.TREG_ACTIVE = self.TREG_ACTIVE
        CC3DType.TREG_ANERGIC = self.TREG_ANERGIC
       
        CC3DType.TCONV_INACTIVE = self.TCONV_INACTIVE
        CC3DType.TCONV_ACTIVE = self.TCONV_ACTIVE
        CC3DType.TCONV_ANERGIC = self.TCONV_ANERGIC
        
    def start(self):   
        # adding options that setup SBML solver integrator - these are optional but useful when encountering integration instabilities              
        # options={'relative':1e-10,'absolute':1e-12}
        # self.setSBMLGlobalOptions(options)
        
        from Cell import APC
        from Cell import TCell
        from Cell import State
        
		
		# Set-up each cell by attaching our own interactions to each one.
		# We use the internal cell dictionary to store information such as concentrations, etc.
        for cell in self.cellList:
			# Get the dictionary of the cell.
			# For more info about the cell dictionary see the CC3D documentation.
			# It's basically map. We always use the same key [CC3DKey.DATA_KEY] to get our information.
            cellDict = self.getDictionaryAttribute(cell)
            
			# Set-up inactive TREGs
            if cell.type == self.TREG_INACTIVE:  
                cellDict[CC3DKey.DATA_KEY] = TCell(cell, TCell.TREG)

			# Set-up active TREGs
            if cell.type == self.TREG_ACTIVE:  
                cellDict[CC3DKey.DATA_KEY] = TCell(cell, TCell.TREG, state=State.ACTIVE)              
            
			# Set-up inactive TCONVs
            elif cell.type == self.TCONV_INACTIVE:
                cellDict[CC3DKey.DATA_KEY] = TCell(cell, TCell.TCONV)
				
				# Attach the recycling SBML file
                #modelFile = 'Simulation/recycling.xml'
				
				# Set-up the initial concentrations
                #initialConditions = {}
                #initialConditions['S1'] = 100
                #initialConditions['S2'] = 1000
                
				# Attach the SBML to the cell
				# The step-size determines how many SBML steps will equal 1 MCS
				# StepSize 1 -> 1 MCS equal to 1 SBML step
				
                #self.addSBMLToCell(_modelFile=modelFile,_modelName='recycling',_cell=cell, _stepSize=1, _initialConditions=initialConditions)
            
			# Set-up active TCONVs
            elif cell.type == self.TCONV_ACTIVE:
                cellDict[CC3DKey.DATA_KEY] = TCell(cell, TCell.TCONV, state=State.ACTIVE)                  
            
			# Set-up APCs
            elif cell.type == self.APC:
                cellDict[CC3DKey.DATA_KEY] = APC(cell)
                
     
    def step(self,mcs):
		# Run the SBML biochemical reaction network. 
        #self.timestepSBML()
        
		for cell in self.cellList:
			############# SBML things are commented out. They require an SBML model named recycling.xml
			############# For more information about implementing BioNetGen and CC3D look at the samples provided by
			############# CC3D that deal with SBMLs in order to simulate CTLA-4 recycling using a biochemical network
			############# I'm currently fixing this part as it isn't as reliable as I wanted it to be.
            #try:
				# Get the current state of the SBML in this cell
                # state = self.getSBMLState(_modelName='recycling',_cell=cell)
				
            #except RuntimeError:
                #pass
            
			# APCs move faster than TCells
			# We try to simulate that here using external potential
            if cell.type == self.APC:
                cell.lambdaVecX = 6
                cell.lambdaVecY = 6
                cell.lambdaVecZ = 6
            else:
                cell.lambdaVecX = 45
                cell.lambdaVecY = 45
                cell.lambdaVecZ = 45
                
            # Otherwise do our usual TCell and APC interaction if they are neighbors (that means they are next to each other)
            for neighbor, commonSurfaceArea in self.getCellNeighborDataList(cell):                
                if neighbor:
                    self.interact(cell, neighbor, mcs)
        
    def interact(self, cell, neighbor, mcs):
		# Get the cell dictionary
        cellDict = self.getDictionaryAttribute(cell)
		# Get the information from that dictionary
        cellInfo = cellDict[CC3DKey.DATA_KEY]
        
		# Get the neighboring cell dictionary
        neighborDict = self.getDictionaryAttribute(neighbor)
		# Get the information from that dictionary
        neighborCellInfo = neighborDict[CC3DKey.DATA_KEY]
        
		# Call our interaction methods depending on the types of cells that are interacting
        if neighbor.type == self.APC:
            cellInfo.interact_with_apc(neighborCellInfo, mcs)
        else:
            cellInfo.interact_with_tcell(neighborCellInfo, mcs)

##########################################################
#	PlotSteppable
#
#	This steppable controls all of the plots displayed in CC3D.
#	When the simulation finishes all plots are saved as pictures in the current working CC3D directory.
#
#	IMPORTANT NOTE:
#		You must let the simulation finish to have the plots be saved.
#		If you STOP the simulation the plots won't be saved.
#		They will only be saved IF the simulation finishes.
##########################################################
class PlotSteppable(SteppableBasePy):
    def __init__(self,_simulator,_frequency=1):
        SteppableBasePy.__init__(self,_simulator,_frequency)    
       
    def start(self):
        # Ligand Losing
        self.pW_lost_ligand = self.addNewPlotWindow(_title='Amount of Peptide-MHC, CD80, CD86 ligand lost',_xAxisTitle='MonteCarlo Step (MCS)',_yAxisTitle='Amount (AU)', _xScaleType='linear',_yScaleType='linear')
        self.pW_lost_ligand.addPlot('Peptide-MHC',_style='Lines',_color='red',_size=2)
        self.pW_lost_ligand.addPlot('CD80',_style='Lines',_color='green',_size=2)
        self.pW_lost_ligand.addPlot('CD86',_style='Lines',_color='purple',_size=2)
        
        # Affinity / Receptor Engagement
        self.pW_affinity = self.addNewPlotWindow(_title='External CTLA-4 and CD28 Receptor Engagement',_xAxisTitle='MonteCarlo Step (MCS)',_yAxisTitle='Number of times bound to ligand', _xScaleType='linear',_yScaleType='linear')
        self.pW_affinity.addPlot('CD28',_style='Lines',_color='blue',_size=2)
        self.pW_affinity.addPlot('External CTLA-4',_style='Lines',_color='red',_size=2)  
        
        # APC Downregulation / Ligands
        self.pW_apc_amount = self.addNewPlotWindow(_title='Peptide-MHC, CD80, and CD86 Downregulation',_xAxisTitle='MonteCarlo Step (MCS)',_yAxisTitle='Amount (AU)', _xScaleType='linear',_yScaleType='linear')
        self.pW_apc_amount.addPlot('Peptide-MHC',_style='Lines',_color='red',_size=2)        
        self.pW_apc_amount.addPlot('CD80',_style='Lines',_color='green',_size=2)
        self.pW_apc_amount.addPlot('CD86',_style='Lines',_color='purple',_size=2)   

        # T-Cell Downregulation / Receptors
        self.pW_tcell_amount = self.addNewPlotWindow(_title='CD28, External CTLA-4, and Internal CTLA-4 Downregulation',_xAxisTitle='MonteCarlo Step (MCS)',_yAxisTitle='Amount (AU)', _xScaleType='linear',_yScaleType='linear')
        self.pW_tcell_amount.addPlot('CD28',_style='Lines',_color='blue',_size=2)
        self.pW_tcell_amount.addPlot('Internal CTLA-4',_style='Lines',_color='red',_size=2)   
        self.pW_tcell_amount.addPlot('External CTLA-4',_style='Lines',_color='green',_size=2) 

        # Total Cell Counts
        self.pW_cell_count = self.addNewPlotWindow(_title='Total Cell Counts',_xAxisTitle='MonteCarlo Step (MCS)',_yAxisTitle='Number of Cells', _xScaleType='linear',_yScaleType='linear')
        self.pW_cell_count.addPlot('APC',_style='Lines',_color='green',_size=2)
        self.pW_cell_count.addPlot('T-Cells',_style='Lines',_color='purple',_size=2)

        # TREG Cell Counts
        self.pW_treg_count = self.addNewPlotWindow(_title='TREG Cells Count',_xAxisTitle='MonteCarlo Step (MCS)',_yAxisTitle='Number of Cells', _xScaleType='linear',_yScaleType='linear')
        self.pW_treg_count.addPlot('Inactive',_style='Lines',_color='green',_size=2)
        self.pW_treg_count.addPlot('Active',_style='Lines',_color='yellow',_size=2)   
        self.pW_treg_count.addPlot('Anergic',_style='Lines',_color='brown',_size=2)

        # TCONV Cell Counts
        self.pW_tconv_count = self.addNewPlotWindow(_title='TCONV Cells Count',_xAxisTitle='MonteCarlo Step (MCS)',_yAxisTitle='Number of Cells', _xScaleType='linear',_yScaleType='linear')
        self.pW_tconv_count.addPlot('Inactive',_style='Lines',_color='magenta',_size=2)
        self.pW_tconv_count.addPlot('Active',_style='Lines',_color='darkblue',_size=2)   
        self.pW_tconv_count.addPlot('Anergic',_style='Lines',_color='cyan',_size=2)

        # Stochastic Occurences
        self.pW_stochastic = self.addNewPlotWindow(_title='Stochastic Processes',_xAxisTitle='MonteCarlo Step (MCS)',_yAxisTitle='Number of Occurences', _xScaleType='linear',_yScaleType='linear')
        self.pW_stochastic.addPlot('Apoptosis',_style='Lines',_color='magenta',_size=2)
        self.pW_stochastic.addPlot('Division',_style='Lines',_color='darkblue',_size=2)   
        self.pW_stochastic.addPlot('Quiescence',_style='Lines',_color='cyan',_size=2)
                       
    def step(self, mcs):
        # Ligand Losing
        self.pW_lost_ligand.addDataPoint('Peptide-MHC', mcs, Data.TOTAL_LOST_PEPTIDEMHC)        
        self.pW_lost_ligand.addDataPoint('CD80', mcs, Data.TOTAL_LOST_CD80)
        self.pW_lost_ligand.addDataPoint('CD86', mcs, Data.TOTAL_LOST_CD86)
        
        # Affinity / Receptor Engagement
        self.pW_affinity.addDataPoint('CD28', mcs, Data.TOTAL_ENGAGED_CD28)
        self.pW_affinity.addDataPoint('External CTLA-4', mcs, Data.TOTAL_ENGAGED_EXTERNAL_CTLA4)
        
        # APC Downregulation / Ligands
        self.pW_apc_amount.addDataPoint('Peptide-MHC', mcs, Data.TOTAL_AMOUNT_PEPTIDEMHC)
        self.pW_apc_amount.addDataPoint('CD80', mcs, Data.TOTAL_AMOUNT_CD80)
        self.pW_apc_amount.addDataPoint('CD86', mcs, Data.TOTAL_AMOUNT_CD86)
        
        # T-Cell Downregulation / Receptors
        self.pW_tcell_amount.addDataPoint('TCR', mcs, Data.TOTAL_AMOUNT_TCR)
        self.pW_tcell_amount.addDataPoint('CD28', mcs, Data.TOTAL_AMOUNT_CD28)
        self.pW_tcell_amount.addDataPoint('External CTLA-4', mcs, Data.TOTAL_AMOUNT_EXTERNAL_CTLA4)
        self.pW_tcell_amount.addDataPoint('Internal CTLA-4', mcs, Data.TOTAL_AMOUNT_INTERNAL_CTLA4) 
        
        # Total Cell Counts
        self.pW_cell_count.addDataPoint('APC', mcs, Data.TOTAL_APC)
        self.pW_cell_count.addDataPoint('T-Cells', mcs, Data.TOTAL_TCELLS)
        
        # TREG Cell Counts
        self.pW_treg_count.addDataPoint('Inactive', mcs, Data.TOTAL_TREG_INACTIVE)
        self.pW_treg_count.addDataPoint('Active', mcs, Data.TOTAL_TREG_ACTIVE)
        self.pW_treg_count.addDataPoint('Anergic', mcs, Data.TOTAL_TREG_ANERGIC)
        
        # TCONV Cell Counts
        self.pW_tconv_count.addDataPoint('Inactive', mcs, Data.TOTAL_TCONV_INACTIVE)
        self.pW_tconv_count.addDataPoint('Active', mcs, Data.TOTAL_TCONV_ACTIVE)
        self.pW_tconv_count.addDataPoint('Anergic', mcs, Data.TOTAL_TCONV_ANERGIC)
        
        # Stochastic Occurences
        self.pW_stochastic.addDataPoint('Apoptosis', mcs, Data.TOTAL_STOCHASTIC_APOPTOSIS)
        self.pW_stochastic.addDataPoint('Division', mcs, Data.TOTAL_STOCHASTIC_DIVISION)
        self.pW_stochastic.addDataPoint('Quiescence', mcs, Data.TOTAL_STOCHASTIC_QUIESCENCE)
            
    def finish(self):
        # Ligand Losing
        self.pW_lost_ligand.savePlotAsData('lost_ligands.dat') 
        self.pW_lost_ligand.savePlotAsPNG('lost_ligands.png', 1000, 1000)
        
        # Affinity / Receptor Engagement
        self.pW_affinity.savePlotAsData('affinities.dat') 
        self.pW_affinity.savePlotAsPNG('affinities.png',1000,1000)
        
        # APC Downregulation / Ligands
        self.pW_apc_amount.savePlotAsData('apc_ligands.dat') 
        self.pW_apc_amount.savePlotAsPNG('apc_ligands.png',1000,1000)
        
        # T-Cell Downregulation / Receptors
        self.pW_tcell_amount.savePlotAsData('tcell_receptors.dat') 
        self.pW_tcell_amount.savePlotAsPNG('tcell_receptors.png',1000,1000)
        
        # Total Cell Counts
        self.pW_cell_count.savePlotAsData('cell_counts.dat') 
        self.pW_cell_count.savePlotAsPNG('cell_counts.png',1000,1000)
        
        # TREG Cell Counts
        self.pW_treg_count.savePlotAsData('treg_cells.dat') 
        self.pW_treg_count.savePlotAsPNG('treg_cells.png',1000,1000)
        
        # TCONV Cell Counts
        self.pW_tconv_count.savePlotAsData('tconv_cells.dat') 
        self.pW_tconv_count.savePlotAsPNG('tconv_cells.png',1000,1000)
        
        # Stochastic Occurences
        self.pW_stochastic.savePlotAsData('stochastic.dat')
        self.pW_stochastic.savePlotAsPNG('stochastic.png',1000,1000)

##########################################################
#	VolumeSteppable
#
#	This steppable allows us to simulate aging of cells.
#	I based this on the examples provided by CC3D for mitosis.
#
##########################################################
class VolumeSteppable(SteppableBasePy):
    def __init__(self,_simulator,_frequency=1):
        SteppableBasePy.__init__(self,_simulator,_frequency)
        
    def start(self):
        for cell in self.cellList:
			# Set our initial age
            cell.targetVolume = Config.INITIAL_AGE
            cell.lambdaVolume = 1.0
            
    def step(self, mcs):
        for cell in self.cellList:
			# Increase our age by the age inside the configuration
            cell.targetVolume += Config.STEP_AGE
            
    def finish(self):
        pass            

##########################################################
#	MitosisSteppable
#
#	This steppable manages what cells do once they reach their decision age.
#	They either undergo apoptosis, division, or quiescence (they don't do anything)
#	You'd be best to not modify this but instead modify the settings inside Config
#	Specifically the Stochastic settings
##########################################################
class MitosisSteppable(MitosisSteppableBase):
    def __init__(self, _simulator, _frequency=1):
        MitosisSteppableBase.__init__(self, _simulator, _frequency)
        
    def step(self, mcs):     
        for cell in self.cellList:
			# Check if the cell has passed our age threshold
            if cell.volume > Config.DECISION_AGE:
				# These are the possible actions the cell might undergo.
                actions_available = ['Apoptosis', 'Division', 'Quiescence']
				# Each action has a corresponding weight or possibility associated with it.
                weights = [Config.PROB_APOPTOSIS, Config.PROB_DIVISION, Config.PROB_QUIESCENCE]
                
				# The numpy library contains a function called choice
				# It allows us to select a choice randomly with a given amount of weights
                from numpy.random import choice
                action = choice(actions_available, p=weights)
                #print 'Stochastically chose ' + str(action) + ' in a cell of type ' + str(cell.type)
                
				# Apoptosis (Cell Death)
                if action == 'Apoptosis':
					# Kill the cell inside CC3D by setting its volume to 0
					# Maybe there's a better way to delete cells.
					# I used one of the CC3D samples for this.
                    cell.targetVolume = 0
					
					# Update our plots
                    Data.TOTAL_STOCHASTIC_APOPTOSIS += 1
                    
					# Decrease the counts depending on what kind of cell it is
                    Data.TOTAL_APC -= 1 if cell.type == self.APC else 0
                    Data.TOTAL_TCELLS -= 1 if cell.type != self.APC else 0
                    
                    Data.TOTAL_TREG_INACTIVE -= 1 if cell.type == self.TREG_INACTIVE else 0
                    Data.TOTAL_TREG_ACTIVE -= 1 if cell.type == self.TREG_ACTIVE else 0
                    Data.TOTAL_TREG_ANERGIC -= 1 if cell.type == self.TREG_ANERGIC else 0
                    
                    Data.TOTAL_TCONV_INACTIVE -= 1 if cell.type == self.TCONV_INACTIVE else 0
                    Data.TOTAL_TCONV_ACTIVE -= 1 if cell.type == self.TCONV_ACTIVE else 0
                    Data.TOTAL_TCONV_ANERGIC -= 1 if cell.type == self.TCONV_ANERGIC else 0
                    
				# Division
				# Not actually implemented
				# Here you would do mitosis as stated in the CC3D mitosis samples
				# We record how many times cells would have divided in a plot
                elif action == 'Division':
                    Data.TOTAL_STOCHASTIC_DIVISION += 1
				
				# Quiescence
				# We don't do anything
				# All we do is increase the count for the plot
                else:
                    Data.TOTAL_STOCHASTIC_QUIESCENCE += 1