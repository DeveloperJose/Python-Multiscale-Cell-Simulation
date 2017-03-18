##########################################################
#	File: MainProgram.py
#	Author: Jose Perez <josegperez@mail.com>
#	Version: Model v5
#
#	This file is part of the CC3D set-up routines
#	
##########################################################

import sys
from os import environ
from os import getcwd
import string

sys.path.append(environ["PYTHON_MODULE_PATH"])

import CompuCellSetup

sim, simthread = CompuCellSetup.getCoreSimulationObjects()
            
CompuCellSetup.initializeSimulationObjects(sim,simthread)
     
steppableRegistry = CompuCellSetup.getSteppableRegistry()
        
from Steppables import MainSteppable
steppableInstance = MainSteppable(sim,_frequency=1)
steppableRegistry.registerSteppable(steppableInstance)

from Steppables import PlotSteppable
steppableInstance = PlotSteppable(sim,_frequency=1)
steppableRegistry.registerSteppable(steppableInstance)

from Steppables import VolumeSteppable
steppableInstance = VolumeSteppable(sim,_frequency=10)
steppableRegistry.registerSteppable(steppableInstance)  

CompuCellSetup.mainLoop(sim,simthread,steppableRegistry)    